function out = GUIcov(x,y,w)
%GUIcov shows the necessary calculations to obtain the covariance in a GUI.
%
%
%<a href="matlab: docsearchFS('GUIcov')">Link to the help function</a>
%
%
%
%
%  Required input arguments:
%
%     x : vector of numeric data or table. Vector or table.
%           Vector containing strictly numerical data.
%           If x is table the second input argument y is not necessary. In
%           this case weighted covariance is computed where the weights are
%           the values inside the contingency table.
%           Data Types - double or table
%
%     y : vector of numeric data. Vector.
%           Vector containing strictly numerical data.
%           This input argument is not requested if previously input
%           argument x is a table.
%           Data Types - double
%
%  Optional input arguments:
%
%    w  : weights. Vector.
%         Vector of the same length of x containing the weights
%         (frequencies) assigned to each observation.
%           Example - 1:10
%           Data Types - double
%
% Output:
%
%    out = detailed output to compute the index. Table. 
%          Table with n+1 rows (where n is the length of x) containing
%          what is shown in the GUI. Last row contains the totals.
%
%
%
%
% See also: GUIvar, GUImad, GUIskewness
%
% References:
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('GUIcov')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%
%{
    % Example of unweighted covariance.
    % x= monthly income of 13 families.
    % y= free time expenditure.
    x=[1330 1225 1225 1400 1575 2050 1750 2240 1225 1730 1470 2730 1380];
    y=[120 60 30 60 90 150 140 210 30 100 30 270 260];
    GUIcov(x,y)
%}

%{
    % Example of calculation of standard deviation (using n as denominator).
    x=[427 492 445 444 476 470];
    GUIstd(x,1)
%}

%{
    % Example 1 of weighted covariance.
    % In this example vectors x y and w are supplied.
    x=[  8     8     9     9];
    y=[6     7     5     7];
    w=[0.4000    0.1000    0.3000    0.2000];
    y=GUIcov(x,y,w);
%}

%{
    %% Example 2 of weighted covariance.
    % In this example first input argument is a table and only this
    % argument is passed.
    N=[0 0.4 0.1
    0.3	 0	0.2];
    colnames={'5' '6'	'7'};
    rownames={'8','9'};
    Ntable=array2table(N,'RowNames',rownames,'VariableNames',colnames);
    out=GUIcov(Ntable);
%}

%% Beginning of code

if istable(x)
    unweighted=false;
    xdouble=str2double(x.Properties.RowNames);
    ydouble=str2double(x.Properties.VariableNames);
    [r,c]=size(x);
    W=table2array(x);
    x=zeros(r*c,1);
    y=x; w=x;
    ij=1;
    for i=1:r
        for j=1:c
            x(ij)=xdouble(i);
            y(ij)=ydouble(j);
            w(ij)=W(i,j);
            ij=ij+1;
        end
    end
    
else
    if nargin<3
        unweighted=true;
    else
        unweighted=false;
    end
    x=x(:);
    y=y(:);
end

lenx=length(x);
seq=(1:lenx)';


if unweighted==true % unweighted standard deviation
    
    sumx=sum(x);
    mx=sumx/lenx;
    sumy=sum(y);
    my=sumy/lenx;
    xmmx=x-mx;
    ymmy=y-my;
    xmmxy=xmmx.*y;
    ymmyx=ymmy.*x;
    numx=sum(xmmxy);
    numy=sum(ymmyx);
    xmmxy=xmmx.*y;
    xmmxymmy=xmmx.*ymmy;
    numxy=sum(xmmxymmy);
    covxy=numxy/lenx;
    header={'i' 'x_i' 'y_i' '(x_i-M_X)' '(y_i-M_Y)' '(x_i - M_X )y_i' '(y_i - M_Y) x_i' '(x_i-M_X)(y_i-M_Y)'};
    
    corpus=[seq, x,y, xmmx, ymmy, xmmxy, ymmyx, xmmxymmy];
    
    footer=[NaN sum(x) sum(y) 0 0, numx, numy, numxy];
    strtitle='Details of covariance $(cov(x,y))$ calculation';
    
else % weighted covariance
    w=w(:);
    n=sum(w);
    sumxw=sum(x.*w);
    mx=sumxw/n;
    sumyw=sum(y.*w);
    my=sumyw/n;
    xmmx=x-mx;
    ymmy=y-my;
    xyw=xmmx.*ymmy.*w;
    
    numxy=sum(xyw);
    covxy=numxy/n;
    header={'i' 'x_i' 'y_i' 'w_i' '(x_i-M_X)' '(y_i-M_Y)' '(x_i-M_X)(y_i-M_Y)w_i'};
    
    corpus=[seq, x,y,w, xmmx, ymmy, xyw];
    
    footer=[NaN sum(x) sum(y) n NaN NaN, numxy];
    strtitle='Details of weighted covariance calculation';
    
end

str=strForSchool(header, corpus, footer);

out=array2table([corpus;footer],'VariableNames',header);

% note that there is a maximum of string size of 1200 characters for the
% LaTeX interpreter
startIndex=regexp(str,'\cr');
a=ceil(length(startIndex)/2);
b=a+1;
toinsert=['\cr ' repmat('&',1,length(header)) ' '];
while length(str) >1150 % 1200
    str=[str(1:startIndex(a)-2) toinsert str(startIndex(b)-1:end)];
    a=a-1;
    b=b-1;
end

fs=14;
dim = [.02 .80 0.1 0.1];
figure('Position',[100 100 1200 600],'Units','normalized');
annotation('textbox',dim,'FitBoxToText','on','String',str,'Interpreter','latex','FontSize',fs);

dim = [.56 .88 0.1 0.1];
if unweighted==true
    strmean=['\boldmath{$M_X$}= $\frac{' num2str(sumx) '}{' num2str(lenx) '}=' num2str(mx)  '\qquad $' ...
        '\boldmath{$M_Y$}= $\frac{' num2str(sumy) '}{' num2str(lenx) '}=' num2str(my) '$'];
else
    strmean=['\boldmath{$M_X$}= $\frac{' num2str(sumxw) '}{' num2str(n) '}=' num2str(mx)  '\qquad $' ...
        '\boldmath{$M_Y$}= $\frac{' num2str(sumyw) '}{' num2str(n) '}=' num2str(my) '$'];
    
end
annotation('textbox',dim,'FitBoxToText','on','String',strmean,'Interpreter','latex','FontSize',fs);

dim = [.02 .9 0.1 0.1];
fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strtitle,'Interpreter','latex','FontSize',fs1);

dim = [.01 .05 0.1 0.1];

% strfin = text at the end of the GUI
if unweighted==true
    strfin=['\boldmath{$cov(x,y)$}=$\frac{\sum_{i=1}^n (x_i-M_X) (y_i-M_Y)}{n}'...
        '=\frac{\sum_{i=1}^n  (x_i-M_X) y_i }{n} = \frac{\sum_{i=1}^n  (y_i-M_Y) x_i }{n}'...
        '= \frac{' num2str(numx) '}{' num2str(lenx) '}=' ...
        num2str(covxy) '$'];
else
    strfin=['\boldmath{$cov(x,y)$}=$\frac{\sum_{i=1}^n (x_i-M_X) (y_i-M_Y)w_i}{\sum_{i=1}^n w_i}'...
        '= \frac{' num2str(numxy) '}{' num2str(n) '}=' ...
        num2str(covxy) '$'];
end

fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strfin,'Interpreter','latex','FontSize',fs1);

end
%FScategory:GUI