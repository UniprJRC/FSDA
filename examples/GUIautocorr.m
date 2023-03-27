function out = GUIautocorr(y, varargin)
%GUIcorr shows the necessary calculations to obtain the autocorrelation in a GUI.
%
%<a href="matlab: docsearchFS('GUIautocorr')">Link to the help function</a>
%
%  Required input arguments:
%
%     y : vector of numeric data or table or timetable.
%           Numeric vector (row or column) or table or timetable.
%           Vector containing strictly numerical data. If x is table or a
%           timetable the first column of the table (timetable) is used to
%           compute the autocorrelation.
%       Data Types - vector of doubles or table or timetable
%
%  Optional input arguments:
%
%    lag  : lag. Positive scalar.
%         Vector of the same length of x containing the weights
%         Positive scalar which contains the required lag to compute the
%         autocorrelation. The default is to use lag 1.
%           Example - 'lag',2
%           Data Types - positive scalar
%
% Output:
%
%    out = detailed output to compute the index. struct.
%          Structure containing the following fields.
%          out.data = table with T+1 rows (where T is the length of y)
%                   containing what is shown in the GUI.
%                   Last row contains the totals.
%          out.corr = scalar containing the autocorrelation coefficient for the chosen lag.
%
%
%
%
% See also: GUIcorr, GUIcov, GUIvar, GUImad, GUIskewness
%
% References:
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('GUIautocorr')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%
%{
    % Details of calculation of autocorrelation of lag 1.
    % Imagine we have the following returns for Microsoft prices:
    % Let us suppose we want to compute the autocovariance at lag 1
    y=[ 5 1 -2 3 -4 6 2 -1 -3 4];
    % If lag is not specified by the default autocorrelation of lag 1 is
    % computed
    out= GUIautocorr(y);

%}

%{
    % Details of calculation of autocorrelation of lag 3.
    % Imagine we have the following returns for Microsoft prices:
    % Let us suppose we want to compute the autocovariance at lag 1
    y=[ 5 1 -2 3 -4 6 2 -1 -3 4];
    out= GUIautocorr(y,'lag',3);
%}


%{
    %% Example of use of  GUIautocorr when input is a timetable
    % Assume y(1) is first day of march 2023 and the series is daily.
    y=[ 5 1 -2 3 -4 6 2 -1 -3 4];
    DateStrings = '2023-03-01';
    t = datetime(DateStrings,'InputFormat','yyyy-MM-dd');
    % Create timetable yt
    yt=array2timetable(y','StartTime',t ,'TimeStep',days);
    % Compute autocorrelation of lag 5
    out= GUIautocorr(yt,'lag',5);
%}



%% Beginning of code
lag=1;
if nargin > 1
    options=struct('lag',lag);

    UserOptions=varargin(1:2:length(varargin));
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:GUIautocorr:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    lag = options.lag;
end


if istable(y) || istimetable(y)
    y=y{:,1};
else
    % Transform y into a column vector
    y=y(:);
end

leny=length(y);
if lag <0 || lag> leny-1
    error('FSDA:GUIautocorr:WrongInputOpt','lag must be an interval in 0, 1, ..., T-1.');
elseif abs(lag-round(lag))>0
    error('FSDA:GUIautocorr:WrongInputOpt','Invalid lag type: must be a real finite integer scalar');
else
end


lagstr=num2str(lag);
ytminuslag=circshift(y,lag);
ytminuslag(1:lag)=NaN;





seq=(1:leny)';

sumy=sum(y);
my=sumy/leny;

ytmy=y-my;
ytminuslagmy=ytminuslag-my;
ytmy2=(y-my).^2;

xmmxymmy=ytmy.*ytminuslagmy;
numcorr=sum(xmmxymmy(lag+1:end));
header={'t' 'y_t' ['y_{t-' lagstr '}'] '(y_t-M_y)' ['(y_{t-' lagstr '}-M_y)'] ...
    ['(y_t-M_y)(y_{t-' lagstr '}-M_y)'] '(y_t - M_y )^2'};

corpus=[seq, y, ytminuslag ytmy, xmmxymmy, ytminuslagmy, ytmy2];
devyt=sum(ytmy2);
ryylagged=numcorr/devyt;

footer=[NaN sum(y) NaN 0 NaN numcorr devyt];
strtitle=['Details of autocorrelation $(corr(y_t,y_{t-' lagstr '}))$ calculation'];

str=strForSchool(header, corpus, footer);


% % note that there is a maximum of string size of 1200 characters for the
% % LaTeX interpreter
% startIndex=regexp(str,'\cr');
% a=ceil(length(startIndex)/2);
% b=a+1;
% % toinsert=['\cr ' repmat('&',1,length(header)) ' '];
% toinsert=['\\ ' repmat('&',1,length(header)) ' '];
% 
% while length(str) >1150 % 1200
%     str=[str(1:startIndex(a)-2) toinsert str(startIndex(b)-1:end)];
%     a=a-1;
%     b=b-1;
% end

fs=14;
dim = [.02 .80 0.1 0.1];
figure('Position',[100 100 1200 600],'Units','normalized');
% Make sure that that figure is also visible inside .mlx files
scatter([],[]);
axis('off')
set(gcf,'Visible','on')
annotation('textbox',dim,'FitBoxToText','on','String',str,'Interpreter','latex','FontSize',fs);

dim = [.56 .88 0.1 0.1];
strmean=['\boldmath{$M_y$}= $\frac{' num2str(sumy) '}{' num2str(leny) '}=' num2str(my) '$'];
annotation('textbox',dim,'FitBoxToText','on','String',strmean,'Interpreter','latex','FontSize',fs);

dim = [.02 .9 0.1 0.1];
fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strtitle,'Interpreter','latex','FontSize',fs1);

dim = [.01 .05 0.1 0.1];

Tstr=num2str(leny);
strfin=['\boldmath{$corr(y_t,y_{t-' lagstr '})$}=$\frac{\sum_{t=' num2str(lag+1) '}^{' Tstr '} (y_t-M_y) (y_{t-' lagstr '}-M_y)}{ \sum_{t=1}^{' Tstr '}(y_t-M_y)^2  }'...
    '= \frac{' num2str(numcorr) '}{' num2str(devyt) '}=' ...
    num2str(ryylagged) '$'];

fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strfin,'Interpreter','latex','FontSize',fs1);
out=struct;
out.data=array2table([corpus;footer],'VariableNames',header);
out.corr=ryylagged;

end
%FScategory:GUI