function out = GUIcorr(x,y,w)
%GUIcorr shows the necessary calculations to obtain the correlation in a GUI.
%
%<a href="matlab: docsearchFS('GUIcorr')">Link to the help function</a>
%
%  Required input arguments:
%
%     x : vector of numeric data or table. Vector or table.
%           Vector containing strictly numerical data. If x is table it
%           contains the contingency table associated with the joint
%           probability distribution. In this case, the second input
%           argument y is not necessary. If x is a table, weighted
%           correation is computed where the weights are the values inside
%           the contingency table. 
%       Data Types - vector of doubles or table
%
%     y : vector of numeric data. Vector.
%           Vector containing strictly numerical data.
%           This input argument is not requested if previously input
%           argument x is a table.
%           Data Types - vector of doubles
%
%  Optional input arguments:
%
%    w  : weights. Vector.
%         Vector of the same length of x containing the weights
%         (frequencies) assigned to each observation.
%           Example - 1:10
%           Data Types - vector of doubles
%
% Output:
%
%    out = detailed output to compute the index. struct.
%          Structure containing the following fields.
%          out.data = table with n+1 rows (where n is the length of x)
%                   containing what is shown in the GUI. 
%                   Last row contains the totals.
%          out.corr = scalar containing the correlation coefficient.
%
%
%
%
% See also: GUIcov, GUIvar, GUImad, GUIskewness, GUIkurtosis
%
% References:
% Milioli, M.A., Riani, M., Zani, S. (2019), "Introduzione all'analisi dei dati statistici (Quarta edizione ampliata)". [MRZ]
% Cerioli, A., Milioli, M.A., Riani, M. (2016), "Esercizi di statistica (Quinta edizione)". [CMR]
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('GUIcorr')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%
%{
    % Example of unweighted covariance.
    % The data below are referred to monthly income of 13 families and
    % their corrisponding free time expenditure (See page 223 of [MRZ]).
    % x= monthly income of 13 families.
    % y= free time expenditure.
    x=[1330 1225 1225 1400 1575 2050 1750 2240 1225 1730 1470 2730 1380];
    y=[120 60 30 60 90 150 140 210 30 100 30 270 260];
    GUIcorr(x,y)
%}

%{
    % Example 1 of weighted correlation.
    % In this example vectors x y and w are supplied. (See Covariance from Wikipedia)
    x=[  8     8     9     9];
    y=[6     7     5     7];
    w=[0.4000    0.1000    0.3000    0.2000];
    y=GUIcorr(x,y,w);
%}

%{
    %% Example 2 of weighted correlation.
    % In this example first input argument is a table and only this
    % argument is passed. (See Covariance from Wikipedia)
    N=[0 0.4 0.1
    0.3	 0	0.2];
    colnames={'5' '6'	'7'};
    rownames={'8','9'};
    Ntable=array2table(N,'RowNames',rownames,'VariableNames',colnames);
    out=GUIcorr(Ntable);
%}

%{
    % Another example of weighted correlation.
    % In this example first input argument is a table and only this
    % argument is passed. (See Correlation from Wikipedia)
    N=[0 1/3 0
    1/3	 0 1/3];
    colnames={'-1' '0' '1'};
    rownames={'0','1'};
    Ntable=array2table(N,'RowNames',rownames,'VariableNames',colnames);
    out=GUIcorr(Ntable);
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
    xmmx2=(x-mx).^2;
    ymmy2=(y-my).^2;
    
    xmmxymmy=xmmx.*ymmy;
    numxy=sum(xmmxymmy);
    header={'i' 'x_i' 'y_i' '(x_i-M_X)' '(y_i-M_Y)' '(x_i - M_X )^2' '(y_i - M_Y)^2' '(x_i-M_X)(y_i-M_Y)'};
    
    corpus=[seq, x,y, xmmx, ymmy, xmmx2, ymmy2, xmmxymmy];
    sumxmmx2=sum(xmmx2); 
    sumymmy2=sum(ymmy2);
        rxy=numxy/sqrt(sumxmmx2*sumymmy2);

    footer=[NaN sum(x) sum(y) 0 0  sumxmmx2 sumymmy2 numxy];
    strtitle='Details of correlation $(corr(x,y))$ calculation';
    
else % weighted correlation
    w=w(:);
    n=sum(w);
    sumxw=sum(x.*w);
    mx=sumxw/n;
    sumyw=sum(y.*w);
    my=sumyw/n;
    xmmx=x-mx;
    ymmy=y-my;
    xyw=xmmx.*ymmy.*w;
    xmmx2=(x-mx).^2.*w;
    ymmy2=(y-my).^2.*w;
    sumxmmx2=sum(xmmx2);
    sumymmy2=sum(ymmy2);

    numxy=sum(xyw);
    rxy=numxy/sqrt(sumxmmx2*sumymmy2);
    header={'i' 'x_i' 'y_i' 'w_i' '(x_i-M_X)^2 w_i' '(y_i-M_Y)^2 w_i' '(x_i-M_X)(y_i-M_Y)w_i'};
    
    corpus=[seq, x,y,w, xmmx2, ymmy2, xyw];
    
    footer=[NaN sum(x) sum(y) n sumxmmx2 sumymmy2, numxy];
    strtitle='Details of weighted correlation calculation';
    
end

str=strForSchool(header, corpus, footer);


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
% Make sure that that figure is also visible inside .mlx files
scatter([],[]);
axis('off')
set(gcf,'Visible','on')
annotation('textbox',dim,'FitBoxToText','on','String',str,'Interpreter','latex','FontSize',fs);

if unweighted==true
dim = [.56 .88 0.1 0.1];
    strmean=['\boldmath{$M_X$}= $\frac{' num2str(sumx) '}{' num2str(lenx) '}=' num2str(mx)  '\qquad $' ...
        '\boldmath{$M_Y$}= $\frac{' num2str(sumy) '}{' num2str(lenx) '}=' num2str(my) '$'];
else
dim = [.46 .90 0.09 0.09];
    strmean=['\boldmath{$M_X$}=$\frac{ \sum_{i=1}^n x_i w_i}{\sum_{i=1}^n w_i}$   = $\frac{' num2str(sumxw) '}{' num2str(n) '}=' num2str(mx)  '\qquad $' ...
             '\boldmath{$M_Y$}=$\frac{ \sum_{i=1}^n y_i w_i}{\sum_{i=1}^n w_i}$= $\frac{' num2str(sumyw) '}{' num2str(n) '}=' num2str(my) '$'];
    
end
annotation('textbox',dim,'FitBoxToText','on','String',strmean,'Interpreter','latex','FontSize',fs);

dim = [.02 .9 0.1 0.1];
fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strtitle,'Interpreter','latex','FontSize',fs1);

dim = [.01 .05 0.1 0.1];

% strfin = text at the end of the GUI
if unweighted==true
    
    strfin=['\boldmath{$corr(x,y)$}=$\frac{\sum_{i=1}^n (x_i-M_X) (y_i-M_Y)}{ \sqrt{\sum_{i=1}^n  (x_i-M_X)^2   \sum_{i=1}^n  (y_i-M_Y)^2 }}'...
        '= \frac{' num2str(numxy) '}{\sqrt{' num2str(sumxmmx2) '\times' num2str(sumymmy2) '}}=' ...
        num2str(rxy) '$'];
else
    strfin=['\boldmath{$corr(x,y)$}=$\frac{\sum_{i=1}^n (x_i-M_X) (y_i-M_Y)w_i}{ \sqrt{\sum_{i=1}^n  (x_i-M_X)^2w_i   \sum_{i=1}^n  (y_i-M_Y)^2 w_i }}'...
        '= \frac{' num2str(numxy) '}{\sqrt{' num2str(sumxmmx2) '\times' num2str(sumymmy2) '}}=' ...
        num2str(rxy) '$'];
end

fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strfin,'Interpreter','latex','FontSize',fs1);
out=struct;
out.data=array2table([corpus;footer],'VariableNames',header);
out.corr=rxy;

end
%FScategory:GUI