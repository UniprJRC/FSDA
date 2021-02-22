function out = GUIregress(x,y, varargin)
%GUIregress shows the necessary calculations to obtain simple linear regression statistics in a GUI.
%
%
%<a href="matlab: docsearchFS('GUIregress')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%     x : vector of numeric data or table. Vector or table.
%           Vector containing strictly numerical data.
%           If x is table the second input argument y is not necessary. In
%           this case the response if the last column of the table.
%           Data Types - double or table
%
%     y : vector of numeric data.
%           Vector containing strictly numerical data.
%           This input argument is not requested if previously input
%           argument x is a table.
%           Data Types - double
%
%  Optional input arguments:
%
%    intercept :  Indicator for constant term. true (default) | false.
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as the comma-separated pair consisting of
%                 'Intercept' and either true to include or false to remove
%                 the constant term from the regression model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
%     plots    : show regression graphically. Boolean.
%                If plots is true an additional plot which shows the (x,y)
%                data, the fitted regression line and the residuals is
%                shown on the screen. Clicking on the line in the legend it
%                is possible to show/hide these three components.
%           Example - 'plots',true
%           Data Types - boolean
%
%    weights  : weights. Vector.
%         Vector of the same length of x containing the weights
%         (frequencies) assigned to each observation.
%           NOT IMPLEMENTED YET
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
%<a href="matlab: docsearchFS('GUIregress')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%
%{
    %% Example of unweighted regression.
    % x= monthly income of 13 families.
    % y= free time expenditure.
    x=[1330 1225 1225 1400 1575 2050 1750 2240 1225 1730 1470 2730 1380];
    y=[120 60 30 60 90 150 140 210 30 100 30 270 260];
    out=GUIregress(x,y);
%}

%{
    %% Use of option plots.
    % The following data matrix reports, for 6 countries, the tourism revenues
    % (y) recorded in a given year (in billions of dollars) and the number of
    % foreign visitors (x) in the same year (in millions of units).
    x=[60 48 45 30 23 15];
    y=[27.3 25.1 58.4 27.1 17.5 11.9];
    out=GUIregress(x,y,'plots',true);
%}

%{
    %% Use of option intercept.
    % The following data matrix reports, for 6 countries, the tourism revenues
    % (y) recorded in a given year (in billions of dollars) and the number of
    % foreign visitors (x) in the same year (in millions of units).
    x=[60 48 45 30 23 15];
    y=[27.3 25.1 58.4 27.1 17.5 11.9];
    out=GUIregress(x,y,'intercept',false);
%}

%{
    % First argument passed as table.
    % The following data matrix reports, for 6 countries, the tourism revenues
    % (y) recorded in a given year (in billions of dollars) and the number of
    % foreign visitors (x) in the same year (in millions of units).
    x=[60 48 45 30 23 15];
    y=[27.3 25.1 58.4 27.1 17.5 11.9];
    XX=array2table([x' y'],'VariableNames',{'x','y'});
    out=GUIregress(XX);
%}


%% Beginning of code

if istable(x)
    y=x{:,end};
    x=x{:,1};
end

lenx=length(x);
seq=(1:lenx)';

weights='';
intercept=true;
plots=false;
if nargin>2
    options=struct('intercept',intercept,'plots',false,'weights',weights);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:GUIregress:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    plots=options.plots;
    intercept=options.intercept;
    weights=options.weights;
end

if isempty(weights)
    unweighted=true;
else
    unweighted=false;
end

if unweighted==true % unweighted regression
    x=x(:);
    y=y(:);
    sumx=sum(x);
    mx=sumx/lenx;
    sumy=sum(y);
    my=sumy/lenx;
    x2=x.^2;
    sumx2=sum(x2);
    xy=x.*y;
    sumxy=sum(xy);
    % achk=numa/denb;
    if intercept ==true
        numb=lenx*sumxy-sumx*sumy;
        denb=lenx*sumx2-sumx^2;
        b=numb/denb;
        a=my-b*mx;
        numa=sumy*sumx2-sumx*sumxy;
    else
        numb=sumxy;
        denb=sumx2;
        b=numb/denb;
        numa=0;
        a=0;
    end
    yhat=a+b*x;
    e=y-yhat;
    sume=sum(e);
    e2=e.^2;
    deve=sum(e2);
    %xe=x.*e;
    sumyhat=sum(yhat);
    
    if intercept==true
        ymhat2=(yhat-sumyhat/lenx).^2;
        ymmy=(y-my).^2;
        devyhat=sum(ymhat2);
        devy=sum(ymmy);
        R2=devyhat/devy;
    else
        ymhat2=yhat.^2;
        ymmy=y.^2;
        devyhat=sum(ymhat2);
        devy=sum(ymmy);
        R2=devyhat/devy;
    end
    
    header={'i' 'x_i' 'y_i' 'x_i^2' 'x_iy_i' '\hat y_i'};
    corpus=[seq, x,y, x2, xy, yhat];
    footer=[NaN sumx sumy sumx2 sumxy sum(yhat)];
    strtitle='Details of calculations of $a$ (intercept),  $b$ (slope)  and $\hat y$ (fitted values)';
    
    if intercept ==true
        header1={'i' 'y_i' '\hat y_i' 'e_i' 'e_i^2' '(\hat y_i -M_Y)^2' '(y_i -M_Y)^2'}; %  'e_i' 'x_ie_i'};
        corpus1=[seq, y , yhat,  e,  e2, ymhat2  ymmy];
        footer1=[NaN sumy sumyhat sume deve devyhat devy];
        strtitle1='Total ($DEV(y)$), regression ($DEV(\hat y)$) and residual deviance ($DEV(e)$)';
    else
        header1={'i' 'y_i' '\hat y_i' 'e_i' 'e_i^2' '\hat y_i^2' 'y_i^2'};
        corpus1=[seq, y , yhat,  e,  e2, ymhat2  ymmy];
        footer1=[NaN sumy sumyhat sume deve devyhat devy];
        strtitle1='Total, regression and residual sum of squares';
    end
else % weighted regression
    % TODO in a future release
    
end

str=strForSchool(header, corpus, footer);
str1=strForSchool(header1, corpus1, footer1);

out=array2table([corpus;footer],'VariableNames',header);

startIndex=regexp(str,'\cr');
a1=ceil(length(startIndex)/2);
b1=a1+1;
toinsert=['\cr ' repmat('&',1,length(header)) ' '];
while length(str) >1150
    str=[str(1:startIndex(a1)-2) toinsert str(startIndex(b1)-1:end)];
    a1=a1-1;
    b1=b1-1;
end

% note that there is a maximum of string size of 1200 characters for the
% LaTeX interpreter
startIndex=regexp(str1,'\cr');
a1=ceil(length(startIndex)/2);
b1=a1+1;
toinsert=['\cr ' repmat('&',1,length(header1)) ' '];
while length(str1) >1100
    str1=[str1(1:startIndex(a1)-2) toinsert str1(startIndex(b1)-1:end)];
    a1=a1-1;
    b1=b1-1;
end

% main window
fs=14;
dim = [.02 .80 0.1 0.1];
h=figure('Position',[100 100 1200 700],'Units','normalized');
h1=figure('Position',[100 100 800 700],'Units','normalized');

annotation(h,'textbox',dim,'FitBoxToText','on','String',str,'Interpreter','latex','FontSize',fs);
drawnow
annotation(h1,'textbox',dim,'FitBoxToText','on','String',str1,'Interpreter','latex','FontSize',fs);

if b>0
    ps='+';
    psm='-';
else
    ps='-';
    psm='+';
end

% Textbox of means and a
if unweighted==true
    strmean=['\boldmath{$M_X$}= $\frac{' num2str(sumx) '}{' num2str(lenx) '}=' num2str(mx)  '\qquad $' ...
        '\boldmath{$M_Y$}= $\frac{' num2str(sumy) '}{' num2str(lenx) '}=' num2str(my) '$'];
    if intercept==true
%         if mx>0
%             signmx='-';
%         else
%             signmx='+';
%         end
        stramean=['\boldmath{$a$}= $M_Y - b M_X =' num2str(my) psm num2str(abs(b)) '\times' num2str(abs(mx)) '=' num2str(a) '$'];
    else
        stramean='\boldmath{$a$}=0';
    end
else
    strmean=['\boldmath{$M_X$}= $\frac{' num2str(sumxw) '}{' num2str(n) '}=' num2str(mx)  '\qquad $' ...
        '\boldmath{$M_Y$}= $\frac{' num2str(sumyw) '}{' num2str(n) '}=' num2str(my) '$'];
end


if sumx>=0
    signsumx='-';
else
    signsumx='+';
end

stryhat1=['\boldmath{$\hat y_1$}= $a+b x_1 =' num2str(a) ps num2str(abs(b)) '\times' num2str(x(1)) '=' num2str(yhat(1)) '$'];
stryhat2=['\boldmath{$\hat y_2$}= $a+b x_2 =' num2str(a) ps num2str(abs(b)) '\times' num2str(x(2)) '=' num2str(yhat(2)) '$'];
stryhati='...';
stryhatn=['\boldmath{$\hat y_{' num2str(lenx)  '}$} = $ a+b x_{' num2str(lenx) '} =' num2str(a) ps num2str(abs(b)) ' \times ' num2str(x(lenx)) '=' num2str(yhat(lenx)) '$'];

posp=0.47;
dimstrmean = [posp .78 0.1 0.1];
annotation(h,'textbox',dimstrmean,'FitBoxToText','on','String',strmean,'Interpreter','latex','FontSize',fs);
dimamean = [posp .68 0.1 0.1];
annotation(h,'textbox',dimamean,'FitBoxToText','on','String',stramean,'Interpreter','latex','FontSize',fs);

dimstryhat1 = [posp .60 0.1 0.1];
annotation(h,'textbox',dimstryhat1,'FitBoxToText','on','String',stryhat1,'Interpreter','latex','FontSize',fs);
dimstryhat2 = [posp .54 0.1 0.1];
annotation(h,'textbox',dimstryhat2,'FitBoxToText','on','String',stryhat2,'Interpreter','latex','FontSize',fs);
dimstryhati = [posp .45 0.1 0.1];
annotation(h,'textbox',dimstryhati,'FitBoxToText','on','String',stryhati,'Interpreter','latex','FontSize',fs);
dimstryhatn = [posp .36 0.1 0.1];
annotation(h,'textbox',dimstryhatn,'FitBoxToText','on','String',stryhatn,'Interpreter','latex','FontSize',fs);

dim = [.02 .9 0.1 0.1];
fs1=20;
annotation(h,'textbox',dim,'FitBoxToText','on','String',strtitle,'Interpreter','latex','FontSize',fs1);
annotation(h1,'textbox',dim,'FitBoxToText','on','String',strtitle1,'Interpreter','latex','FontSize',fs1);

% strfin = text at the end of the GUI
if unweighted==true
    if intercept==true
        strbcoeff=['\boldmath{$b$}=$ \frac{n\sum_{i=1}^n x_i y_i - \sum_{i=1}^n x_i \sum_{i=1}^n y_i}{n \sum_{i=1}^n x_i^2 - \left( \sum_{i=1}^n x_i \right)^2}'...
            '= \frac{' num2str(lenx) ' \times ' num2str(sumxy) signsumx  num2str(abs(sumx)) ' \times ' num2str(sumy)  '}{' num2str(lenx) ' \times ' num2str(sumx2) signsumx num2str(sumx)  '^2}' ...
            '= \frac{' num2str(numb) '}{' num2str(denb) '}=' ...
            num2str(b) '$'];
        stracoeff=['\boldmath{$a$}=$ \frac{\sum_{i=1}^n y_i \sum_{i=1}^nx_i^2 - \sum_{i=1}^n x_i \sum_{i=1}^n x_i y_i}{n \sum_{i=1}^n x_i^2 - \left( \sum_{i=1}^n x_i \right)^2}'...
            '= \frac{' num2str(sumy) ' \times ' num2str(sumx2) signsumx num2str(abs(sumx)) ' \times ' num2str(sumxy)  '}{' num2str(lenx) ' \times ' num2str(sumx2) signsumx num2str(sumx)  '^2}' ...
            '= \frac{' num2str(numa) '}{' num2str(denb) '}=' ...
            num2str(a) '$'];
        
        strR2=['\boldmath{$R^2$}=$ \frac{DEV(\hat{y})}{DEV(y)} '...
            '= \frac{' num2str(devyhat)  '}{' num2str(devy) '}=' ...
            num2str(R2) '$'];
        strR2bis=['\boldmath{$R^2$}=$ 1-\frac{DEV(e)}{DEV(y)}' ...
            '= 1-\frac{' num2str(deve) '}{' num2str(devy) '}=' ...
            num2str(R2) '$'];
    else
        strbcoeff=['\boldmath{$b$}=$ \frac{\sum_{i=1}^n x_i y_i}{ \sum_{i=1}^n x_i^2 }'...
            '= \frac{' num2str(numb) '}{' num2str(denb) '}=' ...
            num2str(b) '$'];
        stracoeff='\boldmath{$a$}=0';
        strR2=['\boldmath{$R^2$}=$ \frac{ \sum_{i=1}^n \hat y_i^2 }{ \sum_{i=1}^n y_i^2} '...
            '= \frac{' num2str(devyhat)  '}{' num2str(devy) '}=' ...
            num2str(R2) '$'];
        strR2bis=['\boldmath{$R^2$}=$ 1-\frac{\sum_{i=1}^n e_i^2 }{\sum_{i=1}^n y_i^2}' ...
            '= 1-\frac{' num2str(deve) '}{' num2str(devy) '}=' ...
            num2str(R2) '$'];
    end
else
    % TODO
end

fs1=20;
dimbcoeff = [.01 .16 0.1 0.1];
annotation(h,'textbox',dimbcoeff,'FitBoxToText','on','String',strbcoeff,'Interpreter','latex','FontSize',fs1);
dimacoeff = [.01 .02 0.1 0.1];
annotation(h,'textbox',dimacoeff,'FitBoxToText','on','String',stracoeff,'Interpreter','latex','FontSize',fs1);

annotation(h1,'textbox',dimbcoeff,'FitBoxToText','on','String',strR2,'Interpreter','latex','FontSize',fs1);
annotation(h1,'textbox',dimacoeff,'FitBoxToText','on','String',strR2bis,'Interpreter','latex','FontSize',fs1);

% Show plot of (x,y) coordinates, fitted line and residuals
if plots==true
    figure
    plot(x,y,'o','LineWidth',2)
    h=refline(b,a);
    %h=lsline;
    h.Color='r';
    hold('on')
    plot([x x]',[yhat y]','k','DisplayName','Residuals')
    if b>0
        abx=['(' num2str(a) '+' num2str(b) 'x)'];
    else
        abx=['(' num2str(a) num2str(b) 'x)'];
    end
    r2str=['(R2=' num2str(R2) ')'];
    leg={'(x,y) coordinates',['Regression line ' abx] ['Residuals ' r2str]};
    clickableMultiLegend(leg,'Location','best')
end

end
%FScategory:GUI