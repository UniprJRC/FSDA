function out = GUItrimmean(x,percent, varargin)
%GUItrimmean shows the necessary calculations to obtain the trimmed mean in a GUI.
%
%
%<a href="matlab: docsearchFS('GUItrimmean')">Link to the help function</a>
%
%
%  Required input arguments:
%
%     x : vector of numeric data. Vector.
%           Vector containg strictly numrical data
%           Data Types - double
%
%    percent  : trimming percent. Scalar.
%           Percentage of input data to be trimmed specified as a scalar
%           between 0 and 100 (generally between 0 and 50).
%           For example, if x is a vector that has n values, the trimmed
%           mean is the mean of x excluding the highest and lowest m data
%           values, where m =[n*(percent/100)/2].
%           For example if n=11 and percent=30, $n*(percent/100)/2=1.65$.
%           $m=[1.65]=1$ therefore the smallest and largest observations are
%           excluded from the computation.
%           Example - 1:10
%           Data Types - double
%
%
%
% Optional input arguments:
%
%
%
% Output:
%
%    out = detailed output to compute the index. Table. Table with n+1 rows (where n is the length of x) containing the
%          what is shown in the GUI. Last row contains the totals.
%
%
%
%
% See also: GUIvar, GUIstd
%
% References:
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('GUItrimmean')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%
%{
    %% Example of use of trimmed mean with 50 per cent of trimming
    x=[60 30 50 50  80 35000  80 95];
    out=GUItrimmean(x,50);
%}

%{
    % Example of use of trimmed mean with 25 per cent of trimming
    x=[60 30 50 50  80 35000  80 95];
    out=GUItrimmean(x,25);
%}

%% Beginning of code

x=x(:);
n=length(x);
seq=(1:n)';

header={'i' 'x_i' 'x_{(i)}' 'x_{(i)} \; not \; trimmed'};

xord=sort(x);
m=floor(n*(percent/100)/2);
xordused=NaN(n,1);
xordused(m+1:n-m)=xord(m+1:n-m);
den=length((m+1):(n-m));
corpus=[seq, x, xord xordused];
num=sum(xord(m+1:n-m));
footer=[NaN sum(x) sum(x) num];
strtitle=['Details of trimmed mean calculation $\alpha=$' num2str(percent/100) ' $m$=' num2str(m)];

str=strForSchool(header, corpus, footer);

out=array2table([corpus;footer],'VariableNames',header);

fs=14;
dim = [.2 .80 0.1 0.1];
figure('Position',[100 100 1000 600],'Units','normalized');
annotation('textbox',dim,'FitBoxToText','on','String',str,'Interpreter','latex','FontSize',fs);

dim = [.2 .9 0.1 0.1];
fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strtitle,'Interpreter','latex','FontSize',fs1);

dim = [.2 .05 0.1 0.1];

trimean=num/den;
nminusmstr=num2str(n-m);
nminus2m=num2str(n-2*m);
str1=['\frac{\sum_{i='  num2str(m+1) '}^{' nminusmstr '} x_{(i)}}{' nminus2m '}='];

strfin=[' \it M$_\alpha(X)= \frac{\sum_{i=m+1}^{n-m} x_{(i)}}{n-2m}='   str1 '\frac{'  num2str(num) '}{' num2str(den) '}=' num2str(trimean) '$'];

fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strfin,'Interpreter','latex','FontSize',fs1);

end
