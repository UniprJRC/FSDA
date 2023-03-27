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
%     x : vector of numeric data or table. Vector or table or empty.
%           Vector containing strictly numerical data.
%           If x is table the second input argument y is not necessary. In
%           this case the response is the last column of the table.
%           If x is empty it is assumed that regression is against time.
%           Data Types - double or table or []
%
%     y : vector of numeric data.
%           Vector containing strictly numerical data.
%           This input argument is not requested if previous input
%           argument x is a table.
%           Data Types - double
%
%  Optional input arguments:
%
%  inferential :  this input parameter controls whether to show an additional GUI
%                 giving all the details for the calculation of all the
%                 inferential aspects of the linear regression model (i.e.
%                 standard error of coefficients, t-statistics, p-values
%                 of the t-statistics for a two-sided test with the null
%                 hypothesis that the coefficient is zero and confidence
%                 intervals). If inferential is empty [] (default) this
%                 additional GUI is now shown. If inferential is a number
%                 in the interval (0 1), this GUI is shown and the a
%                 confidence level at inferential-per cent level is
%                 computed. For example if inferential is 0.90 a 90 per
%                 cent confidence interval is computed.
%                 Example - 'inferential',0.99
%                 Data Types - boolean
%
%    intercept :  Indicator for constant term. true (default) | false.
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as the comma-separated pair consisting of
%                 'Intercept' and either true to include or false to remove
%                 the constant term from the regression model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
%    timeseries : Flag indicating that data is a time series. false (default) | true.
%                 Example - 'timeseries', true
%                 Data Types - boolean
%
%
%     plots    : show regression graphically. Boolean.
%                If plots is true an additional plot which shows the (x,y)
%                data, the fitted regression line and the residuals is
%                shown on the screen. Clicking on the line in the legend it
%                is possible to show/hide these three components.
%           Example - 'plots',true
%           Data Types - boolean
%
%   interpolant  : interpolation method. 'character'.
%         Possible values are 'linear' (default) for linear interpolation,
%         'exponential' for exponential interpolation and 'power' for power
%         interpolation.
%           Example - 'interpolant','exponential'
%           Data Types - char
%
% Output:
%
%
%         out:   structure which contains the following fields
%
%    out.tabledata = detailed output to compute the index. Table.
%          Table with n+1 rows (where n is the length of x) containing
%          what is shown in the GUI. Last row contains the totals.
%    out.a = intercept of linear/exponential/power regression on the
%           original data. scalar.
%     out.b = slope of the of linear/exponential/power regression on the original data. scalar.
%    out.expa = intercept of the exponential/power regression on the
%               transformed data. scalar. This field is present only if
%               input option interpolant is 'exponential' or 'power'.
%    out.expb =  slope of the exponential regression on the transformed data. scalar.
%                This field is present only if input option interpolant is 'exponential'.
%    out.R2 = R squared of the regression model. scalar.
%            When exponential/power regression is chosen, R squared is always
%            computed on transformed data.
%    out.ta  = t statistic for intercept. Scalar.
%    out.tb  = t statistic for slope. Scalar.
%    out.confinta = confidence interval for alpha parameter. A vector of
%            length 2 containing the extremes of confidence interval for the
%            intercept. The default confidence level is 0.95 or it
%            is controlled by optional input argument inferential.
%    out.confintb = confidence interval for alpha parameter. A vector of
%            length 2 containing the extremes of confidence interval for the
%            slope. The default confidence level is 0.95 or it
%            is controlled by optional input argument inferential.
%
% See also: GUIvar, GUImad, GUIskewness
%
% References:
% Milioli, M.A., Riani, M., Zani, S. (2019), "Introduzione all'analisi dei dati statistici (Quarta edizione ampliata)". [MRZ]
% Cerioli, A., Milioli, M.A., Riani, M. (2016), "Esercizi di statistica (Quinta edizione)". [CMR]
%
% Copyright 2008-2023.
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
    %% Calculation of unweighted regression.
    % In this example we know the monthly income of 13 families and we
    % estimate the correlation with the free time expenditure. (See page 223 of [MRZ]).
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
    % foreign visitors (x) in the same year (in millions of units). (See page 101 of [CMR])
    x=[60 48 45 30 23 15];
    y=[27.3 20.1 25.4 27.1 17.5 11.9];
    out=GUIregress(x,y,'plots',true);
%}

%{
    %% Use of option intercept.
    % The following data matrix reports, for 6 countries, the tourism revenues
    % (y) recorded in a given year (in billions of dollars) and the number of
    % foreign visitors (x) in the same year (in millions of units). (See page 101 of [CMR])
    x=[60 48 45 30 23 15];
    y=[27.3 25.1 58.4 27.1 17.5 11.9];
    out=GUIregress(x,y,'intercept',false,'plots',false);
%}

%{
    % First argument passed as table.
    % The following data matrix reports, for 6 countries, the tourism revenues
    % (y) recorded in a given year (in billions of dollars) and the number of
    % foreign visitors (x) in the same year (in millions of units). (See page 101 of [CMR])
    x=[60 48 45 30 23 15];
    y=[27.3 25.1 58.4 27.1 17.5 11.9];
    XX=array2table([x' y'],'VariableNames',{'x','y'});
    out=GUIregress(XX);
%}

%{
    %% Example of exponential interpolation.
    % The values ​​of a company's production, in millions of euros were as follows: (See page 116 of [CMR])  
    x=1:7;
    y=[50 52 56 59 65 70 76]; 
    % Analyze the trend of the company's production using an exponential fit.
    out=GUIregress([],y,'interpolant','exponential','plots',true,'timeseries', true);
%}

%{
    %% Example of power interpolation.
    % The values ​​of a company's production, in millions of euros were as follows: (See page 116 of [CMR])   
    x=1:7;
    y=[50 52 56 59 65 70 76];
    % Analyze the trend of the company's production using a power fit.
    out=GUIregress(x,y,'interpolant','power','plots',true,'timeseries', true);
%}


%{
    %% Example of linear, exponential and power interpolation .
    % Time series of the value ​​of a commodity, in euros were as follows: (See page 269 of [MRZ])   
    y=[50496 52396 55058 56550 56275 58138 59485 61188 63989 66505];
    x=1:10;
    % Analyze the trend of the company's production using a linear fit.
    out=GUIregress(x,y,'interpolant','linear','plots',true, 'timeseries', true);

    % Analyze the trend of the company's production using an exponential fit.
    out=GUIregress(x,y,'interpolant','exponential','plots',true,'timeseries',true);

    % Analyze the trend of the company's production using an power fit.
    out=GUIregress(x,y,'interpolant','power','plots',true,'timeseries',true)
    

%}

%{
    %% Example of linear and power interpolation .
    close all
    % Time series....,  (See ex 4.26 of [CMR])   
    xb=1:8;
    xa=[1 3 5 7 9 11 13 15];
    y=[325 327 329 332 335 338 340 343];

    % Analyze the trend of the company's production using a linear fit.
    out=GUIregress(xa,y,'interpolant','linear','plots',true, 'timeseries', true);

    % Analyze the trend of the company's production using an exponential fit.
    out=GUIregress(xb,y,'interpolant','exponential','plots',true, 'timeseries', true);

     % Analyze the trend of the company's production using a power fit.
    out=GUIregress(xa,y,'interpolant','power','plots',true, 'timeseries', true);
%}

%{
    %% Use of option inferential.
    % In a survey on pollution, we want to verify whether the content of a
    % certain pollutant in the air, expressed in micrograms per cubic meter
    % (Y) is linked to the number of manufacturing companies with more than
    % 20 employees (X). The results obtained in some cities are given
    % below:
    x=[91 453 254 412 334 428 341 125];
    y=[13 12 17 56 29 35 49 27];
    out=GUIregress(x,y,'inferential',0.90,'plots',false);
    % There is not enough evidence to state that is not 0 in the population
    % and the estimated relationship between the two variables is
    % significant. The value of the linear index of determination
    % (R2=0.159) shows a very poor fit of the regression line. The
    % estimated model is therefore of no use for interpreting the causes of
    % pollution.
%}

%% Beginning of code

if isempty(x)
    x=1:length(y);
    timeseries=true;

elseif istable(x)
    y=x{:,end};
    x=x{:,1};
    timeseries=false;
else
    timeseries=false;
end

lenx=length(x);
seq=(1:lenx)';

% weights='';
% %    weights  : weights. Vector.
% %         Vector of the same length of x containing the weights
% %         (frequencies) assigned to each observation.
% %           NOT IMPLEMENTED YET
% %           Example - 1:10
% %           Data Types - double

intercept=true;
plots=false;
interpolant='linear';
inferential=[];
if nargin>2
    options=struct('intercept',intercept,'plots',false,...
        'interpolant',interpolant, 'timeseries', false,'inferential',inferential);

    % 'weights',weights,...

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
    % weights=options.weights;
    interpolant=options.interpolant;
    timeseries=options.timeseries;
    inferential=options.inferential;
end

weights=[];
if isempty(weights)
    unweighted=true;
else
    unweighted=false;
end

if unweighted==true % unweighted regression
    x=x(:);
    y=y(:);

    if strcmp(interpolant,'exponential')
        yori=y;
        y=log(y);
        interpType=2;
    elseif strcmp(interpolant,'power')
        xori=x;
        yori=y;
        x=log(x);
        y=log(y);
        interpType=3;
    else
        interpType=1;
    end

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

    if timeseries==true
        if  interpType==1
            header={'t' 'y_t' 't^2' 'ty_t' '\hat y_t'};
            corpus=[x,y, x2, xy, yhat];
            footer=[sumx sumy sumx2 sumxy sum(yhat)];
        elseif interpType==2
            yhat=exp(yhat);
            header={'t' 'log(y_t)' 't^2' 'tlog(y_t)' 'y_t' '\hat y_t'};
            corpus=[x,y, x2, xy, yori, yhat];
            footer=[sumx sumy sumx2 sumxy sum(yori) sum(yhat)];
        elseif interpType==3
            yhat=exp(a)*xori.^b;
            header={'log(t)' 'log(y_t)' 'log(t)^2' 'log(t)log(y_t)' 'y_t' '\hat y_t'};
            corpus=[x,y,           x2, xy, yori, yhat];
            footer=[sumx sumy sumx2 sumxy sum(yori) sum(yhat)];
        else
        end

    else
        header={'i' 'x_i' 'y_i' 'x_i^2' 'x_iy_i' '\hat y_i'};
        corpus=[seq, x,y, x2, xy, yhat];
        footer=[NaN sumx sumy sumx2 sumxy sum(yhat)];
    end
    strtitle='Details of calculations of $a$ (intercept),  $b$ (slope)  and $\hat y$ (fitted values)';

    if intercept==true
        if timeseries==true
            if interpType==1 % linear
                header1={'t' 'y_t' '\hat y_t' 'e_t' 'e_t^2' '(\hat y_t -M_Y)^2' '(y_t -M_Y)^2'};
                corpus1=[seq, y , yhat,  e,  e2, ymhat2  ymmy];
                footer1=[NaN sumy sumyhat sume deve devyhat devy];
            elseif interpType==2 % exponential
                header1={'t' 'log y_t' 'log \hat y_t' 'e_t' 'e_t^2' '( log \hat y_t -M_Y)^2' '( log y_t -M_Y)^2'};
                corpus1=[seq, y , log(yhat),  e,  e2, ymhat2  ymmy];
                footer1=[NaN sumy sum(log(yhat)) sume deve devyhat devy];
            elseif interpType==3 % power
                header1={'t' 'log t' 'log y_t' 'log \hat y_t' 'e_t' 'e_t^2' '( log \hat y_t -M_Y)^2' '( log y_t -M_Y)^2'};
                corpus1=[seq, x, y , log(yhat),  e,  e2, ymhat2  ymmy];
                footer1=[NaN sumx sumy sum(log(yhat)) sume deve devyhat devy];
            else
            end
        else
            header1={'i' 'y_i' '\hat y_i' 'e_i' 'e_i^2' '(\hat y_i -M_Y)^2' '(y_i -M_Y)^2'};
            corpus1=[seq, y , yhat,  e,  e2, ymhat2  ymmy];
            footer1=[NaN sumy sumyhat sume deve devyhat devy];
        end
        strtitle1='Total ($DEV(y)$), regression ($DEV(\hat y)$) and residual deviance ($DEV(e)$)';
    else % no intercept
        if timeseries==true
            if interpType==1 % linear
                header1={'t' 'y_t' '\hat y_t' 'e_t' 'e_t^2' '\hat y_t^2' 'y_t^2'};
                corpus1=[seq, y , yhat,  e,  e2, ymhat2  ymmy];
                footer1=[NaN sumy sumyhat sume deve devyhat devy];
            elseif interpType==2 % exponential
                header1={'t' 'log(y_t)' '\hat y_t' 'e_t' 'e_t^2' '\hat y_t^2' 'y_t^2'};
                corpus1=[seq, y , yhat,  e,  e2, ymhat2  ymmy];
                footer1=[NaN sumy sumyhat sume deve devyhat devy];
            elseif interpType==3 % power
                header1={'t' 'log t' 'log(y_t)' '\hat y_t' 'e_t' 'e_t^2' '\hat y_t^2' 'y_t^2'};
                corpus1=[seq, x, y , yhat,  e,  e2, ymhat2  ymmy];
                footer1=[NaN sumx sumy sumyhat sume deve devyhat devy];
            end
        else
            header1={'i' 'y_i' '\hat y_i' 'e_i' 'e_i^2' '\hat y_i^2' 'y_i^2'};
            corpus1=[seq, y , yhat,  e,  e2, ymhat2  ymmy];
            footer1=[NaN sumy sumyhat sume deve devyhat devy];
        end

        strtitle1='Total, regression and residual sum of squares';
    end
else % weighted regression
    % TODO in a future release

end

str=strForSchool(header, corpus, footer);
str1=strForSchool(header1, corpus1, footer1);


out = struct;

out.tabledata=array2table([corpus;footer],'VariableNames',header);
if interpType == 1
    out.a=a;
    out.b=b;
    out.R2=R2;
elseif interpType == 2
    out.expa=a;
    out.expb=b;
    out.a=exp(a);
    out.b=exp(b);
    out.R2=R2;
else
    out.expa=a;
    out.a=exp(a);
    out.b=b;
    out.R2=R2;
end

% startIndex=regexp(str,'\cr');
% a1=ceil(length(startIndex)/2);
% b1=a1+1;
% toinsert=['\cr ' repmat('&',1,length(header)) ' '];
% while length(str) >1150
%     str=[str(1:startIndex(a1)-2) toinsert str(startIndex(b1)-1:end)];
%     a1=a1-1;
%     b1=b1-1;
% end
% 
% % note that there is a maximum of string size of 1200 characters for the
% % LaTeX interpreter
% startIndex=regexp(str1,'\cr');
% a1=ceil(length(startIndex)/2);
% b1=a1+1;
% toinsert=['\cr ' repmat('&',1,length(header1)) ' '];
% while length(str1) >1100
%     str1=[str1(1:startIndex(a1)-2) toinsert str1(startIndex(b1)-1:end)];
%     a1=a1-1;
%     b1=b1-1;
% end

% main figure
h=figure('Position',[100 100 1200 700],'Units','normalized','Name','Parameter estimation and fitted values');
h1=figure('Position',[100 100 900 700],'Units','normalized','Name','Deviance decomposition and R2');

% Textboxes for statistics calculation
fs=14;
dim = [.02 .830 0.1 0.1];
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
if timeseries==true
    if interpType==1
        stryhat1=['\boldmath{$\hat y_1$}= $a+b  =' num2str(a) ps num2str(abs(b)) '\times' num2str(x(1)) '=' num2str(yhat(1)) '$'];
        stryhat2=['\boldmath{$\hat y_2$}= $a+b \times 2 =' num2str(a) ps num2str(abs(b)) '\times' num2str(x(2)) '=' num2str(yhat(2)) '$'];
        stryhatn=['\boldmath{$\hat y_{' num2str(lenx)  '}$} = $ a+b \times'  num2str(lenx) ' =' num2str(a) ps num2str(abs(b)) ' \times ' num2str(x(lenx)) '=' num2str(yhat(lenx)) '$'];

    elseif interpType==2 % exponential function
        stryhat1=['\boldmath{$\hat y_1$}= $\exp(a+b  )= \exp(' num2str(a) ps num2str(abs(b)) '\times' num2str(x(1)) ')=' num2str(yhat(1)) '$'];
        stryhat2=['\boldmath{$\hat y_2$}= $\exp(a+b \times 2) = \exp(' num2str(a) ps num2str(abs(b)) '\times' num2str(x(2)) '=' num2str(yhat(2)) '$'];
        stryhatn=['\boldmath{$\hat y_{' num2str(lenx)  '}$} = $ \exp(a+b \times'  num2str(lenx) ') = \exp(' num2str(a) ps num2str(abs(b)) ' \times ' num2str(x(lenx)) ')=' num2str(yhat(lenx)) '$'];
    elseif interpType==3   % power function
        stryhat1=['\boldmath{$\hat y_1$}= $\exp(a)  = \exp(' num2str(a) ') =' num2str(yhat(1)) '$'];
        stryhat2=['\boldmath{$\hat y_2$}= $\exp(a) \times 2^b  = \exp(' num2str(a) ') \times 2^{' num2str(b)  '}=' num2str(yhat(2)) '$'];
        stryhatn=['\boldmath{$\hat y_{' num2str(lenx)  '}$} = $ \exp(a)'  num2str(xori(lenx)) '^{ '  num2str(b) '}=' num2str(yhat(lenx)) '$'];
    end
else
    stryhat1=['\boldmath{$\hat y_1$}= $a+b x_1 =' num2str(a) ps num2str(abs(b)) '\times' num2str(x(1)) '=' num2str(yhat(1)) '$'];
    stryhat2=['\boldmath{$\hat y_2$}= $a+b x_2 =' num2str(a) ps num2str(abs(b)) '\times' num2str(x(2)) '=' num2str(yhat(2)) '$'];
    stryhatn=['\boldmath{$\hat y_{' num2str(lenx)  '}$} = $ a+b x_{' num2str(lenx) '} =' num2str(a) ps num2str(abs(b)) ' \times ' num2str(x(lenx)) '=' num2str(yhat(lenx)) '$'];
end
stryhati='...';

posp=0.49;
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
        if timeseries==true
            if interpType==1  % linear
                strbcoeff=['\boldmath{$b$}=$ \frac{T\sum_{t=1}^T t y_i - \sum_{t=1}^T t \sum_{t=1}^T y_i}{T \sum_{t=1}^T t^2 - \left( \sum_{t=1}^T t \right)^2}'...
                    '= \frac{' num2str(lenx) ' \times ' num2str(sumxy) signsumx  num2str(abs(sumx)) ' \times ' num2str(sumy)  '}{' num2str(lenx) ' \times ' num2str(sumx2) signsumx num2str(sumx)  '^2}' ...
                    '= \frac{' num2str(numb) '}{' num2str(denb) '}=' ...
                    num2str(b) '$'];
                stracoeff=['\boldmath{$a$}=$ \frac{\sum_{t=1}^n y_i \sum_{i=1}^n t^2 - \sum_{t=1}^n t \sum_{i=1}^n t y_i}{T \sum_{t=1}^n t^2 - \left( \sum_{t=1}^T t \right)^2}'...
                    '= \frac{' num2str(sumy) ' \times ' num2str(sumx2) signsumx num2str(abs(sumx)) ' \times ' num2str(sumxy)  '}{' num2str(lenx) ' \times ' num2str(sumx2) signsumx num2str(sumx)  '^2}' ...
                    '= \frac{' num2str(numa) '}{' num2str(denb) '}=' ...
                    num2str(a) '$'];
            elseif interpType==2 % exponential
                strbcoeff=['\boldmath{$b$}=$ \frac{T\sum_{t=1}^T t log y_t - \sum_{t=1}^T t \sum_{t=1}^T log y_t}{T \sum_{t=1}^T t^2 - \left( \sum_{t=1}^T t \right)^2}'...
                    '= \frac{' num2str(lenx) ' \times ' num2str(sumxy) signsumx  num2str(abs(sumx)) ' \times ' num2str(sumy)  '}{' num2str(lenx) ' \times ' num2str(sumx2) signsumx num2str(sumx)  '^2}' ...
                    '= \frac{' num2str(numb) '}{' num2str(denb) '}=' ...
                    num2str(b) '$'];
                stracoeff=['\boldmath{$a$}=$ \frac{\sum_{t=1}^n log y_t \sum_{i=1}^n t^2 - \sum_{t=1}^n t \sum_{i=1}^n t log y_t}{T \sum_{t=1}^n t^2 - \left( \sum_{t=1}^T t \right)^2}'...
                    '= \frac{' num2str(sumy) ' \times ' num2str(sumx2) signsumx num2str(abs(sumx)) ' \times ' num2str(sumxy)  '}{' num2str(lenx) ' \times ' num2str(sumx2) signsumx num2str(sumx)  '^2}' ...
                    '= \frac{' num2str(numa) '}{' num2str(denb) '}=' ...
                    num2str(a) '$'];

            elseif  interpType==3  % power
                strbcoeff=['\boldmath{$b$}=$ \frac{T\sum_{t=1}^T log(t) log y_t - \sum_{t=1}^T log(t) \sum_{t=1}^T log y_t}{T \sum_{t=1}^T log(t)^2 - \left( \sum_{t=1}^T log(t) \right)^2}'...
                    '= \frac{' num2str(lenx) ' \times ' num2str(sumxy) signsumx  num2str(abs(sumx)) ' \times ' num2str(sumy)  '}{' num2str(lenx) ' \times ' num2str(sumx2) signsumx num2str(sumx)  '^2}' ...
                    '= \frac{' num2str(numb) '}{' num2str(denb) '}=' ...
                    num2str(b) '$'];
                stracoeff=['\boldmath{$a$}=$ \frac{\sum_{t=1}^T log y_t \sum_{i=1}^T log(t)^2 - \sum_{t=1}^T log(t) \sum_{i=1}^T log(t) log y_t}{T \sum_{t=1}^n log(t)^2 - \left( \sum_{t=1}^T log(t) \right)^2}'...
                    '= \frac{' num2str(sumy) ' \times ' num2str(sumx2) signsumx num2str(abs(sumx)) ' \times ' num2str(sumxy)  '}{' num2str(lenx) ' \times ' num2str(sumx2) signsumx num2str(sumx)  '^2}' ...
                    '= \frac{' num2str(numa) '}{' num2str(denb) '}=' ...
                    num2str(a) '$'];


            end
        else
            strbcoeff=['\boldmath{$b$}=$ \frac{n\sum_{i=1}^n x_i y_i - \sum_{i=1}^n x_i \sum_{i=1}^n y_i}{n \sum_{i=1}^n x_i^2 - \left( \sum_{i=1}^n x_i \right)^2}'...
                '= \frac{' num2str(lenx) ' \times ' num2str(sumxy) signsumx  num2str(abs(sumx)) ' \times ' num2str(sumy)  '}{' num2str(lenx) ' \times ' num2str(sumx2) signsumx num2str(sumx)  '^2}' ...
                '= \frac{' num2str(numb) '}{' num2str(denb) '}=' ...
                num2str(b) '$'];
            stracoeff=['\boldmath{$a$}=$ \frac{\sum_{i=1}^n y_i \sum_{i=1}^nx_i^2 - \sum_{i=1}^n x_i \sum_{i=1}^n x_i y_i}{n \sum_{i=1}^n x_i^2 - \left( \sum_{i=1}^n x_i \right)^2}'...
                '= \frac{' num2str(sumy) ' \times ' num2str(sumx2) signsumx num2str(abs(sumx)) ' \times ' num2str(sumxy)  '}{' num2str(lenx) ' \times ' num2str(sumx2) signsumx num2str(sumx)  '^2}' ...
                '= \frac{' num2str(numa) '}{' num2str(denb) '}=' ...
                num2str(a) '$'];
        end
        strR2=['\boldmath{$R^2$}=$ \frac{DEV(\hat{y})}{DEV(y)} '...
            '= \frac{' num2str(devyhat)  '}{' num2str(devy) '}=' ...
            num2str(R2) '$'];
        strR2bis=['\boldmath{$R^2$}=$ 1-\frac{DEV(e)}{DEV(y)}' ...
            '= 1-\frac{' num2str(deve) '}{' num2str(devy) '}=' ...
            num2str(R2) '$'];
    else % intercept = false
        if timeseries==true
            if interpType==1
                strbcoeff=['\boldmath{$b$}=$ \frac{\sum_{t=1}^T t y_t}{ \sum_{t=1}^T t^2 }'...
                    '= \frac{' num2str(numb) '}{' num2str(denb) '}=' ...
                    num2str(b) '$'];
            elseif interpType==2
                strbcoeff=['\boldmath{$b$}=$ \frac{\sum_{t=1}^T t log(y_t)}{ \sum_{t=1}^T t^2 }'...
                    '= \frac{' num2str(numb) '}{' num2str(denb) '}=' ...
                    num2str(b) '$'];
            elseif interpType==3
                strbcoeff=['\boldmath{$b$}=$ \frac{\sum_{t=1}^T log(t) log(y_t)}{ \sum_{t=1}^T log(t)^2 }'...
                    '= \frac{' num2str(numb) '}{' num2str(denb) '}=' ...
                    num2str(b) '$'];
            end
        else
            strbcoeff=['\boldmath{$b$}=$ \frac{\sum_{i=1}^n x_i y_i}{ \sum_{i=1}^n x_i^2 }'...
                '= \frac{' num2str(numb) '}{' num2str(denb) '}=' ...
                num2str(b) '$'];
        end
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


% Textboxes for b= and a=
fs1=20;
dimbcoeff = [.01 .16 0.1 0.1];
dimacoeff = [.01 .02 0.1 0.1];
annotation(h,'textbox',dimacoeff,'FitBoxToText','on','String',stracoeff,'Interpreter','latex','FontSize',fs1);

annotation(h,'textbox',dimbcoeff,'FitBoxToText','on','String',strbcoeff,'Interpreter','latex','FontSize',fs1);

% Textboxes for R^2
dimR2coeff = [.01 .16 0.1 0.1];
dimR2biscoeff = [.01 .02 0.1 0.1];
annotation(h1,'textbox',dimR2coeff,'FitBoxToText','on','String',strR2,'Interpreter','latex','FontSize',fs1);
annotation(h1,'textbox',dimR2biscoeff,'FitBoxToText','on','String',strR2bis,'Interpreter','latex','FontSize',fs1);

xmx2=(x-mx).^2;
sumxmx2=sum(xmx2);
s=sqrt(deve/(lenx-2));
stderra=s*sqrt(1/lenx+mx^2/sumxmx2);
stderrb=s*sqrt(1/sumxmx2);
ta=a/stderra;
tb=b/stderrb;
pvala=2*tcdf(abs(ta),lenx-2,'upper');
pvalb=2*tcdf(abs(tb),lenx-2,'upper');
alpha=1-inferential;
tinvinf=tinv(alpha/2,lenx-2);
tinvsup=tinv(1-alpha/2,lenx-2);
tinfa=a+tinvinf*stderra;
tsupa=a+tinvsup*stderra;
tinfb=b+tinvinf*stderrb;
tsupb=b+tinvsup*stderrb;

if intercept == true
    out.ta=ta;
else
    out.ta=[];
end
out.tb=tb;
if intercept == true
    out.confinta=[tinfa tsupa];
else
   out.confinta=[];
end
out.confintb=[tinfb tsupb];


if ~isempty(inferential)
    h2=figure('Position',[100 100 1200 700],'Units','normalized','Name','Inferential aspects');

    header2={'i' '(x_i -M_X)^2'};
    corpus2=[seq, xmx2];
    footer2=[NaN sumxmx2];
    str2=strForSchool(header2, corpus2, footer2);
    annotation(h2,'textbox',dim,'FitBoxToText','on','String',str2,'Interpreter','latex','FontSize',fs);




    strs=['\boldmath{$s= \sqrt{\frac{DEV(E)}{n-2}} $}= $\sqrt{\frac{' num2str(deve) '}{' num2str(lenx-2) '}}=' num2str(s)  '$'];

    strstderrb=['\boldmath{$s_b= s \sqrt{\frac{1}{\sum_{i=1}^n (x_i -M_X)^2}} = '  num2str(s)   '\sqrt{\frac{' num2str(1) '}{' num2str(sumxmx2) '}}= $ $' num2str(stderrb) '$}' ];

    if intercept==true
        strstderra=['\boldmath{$s_a= s \sqrt{\frac{1}{n} +\frac{\overline{x}^2}{\sum_{i=1}^n (x_i -M_X)^2}} $= $'  num2str(s)   '\sqrt{\frac{' num2str(1) '}{' num2str(lenx) '}+\frac{' num2str(mx) '^2}{' num2str(sumxmx2) '}}= $ $' ...
            num2str(stderra) '$}' ];
    else
        strstderra='\boldmath{$NA$}';
    end


    posp=0.39;

    % Root mean squared error
    dim = [posp .9 0.1 0.1];
    fs1=20;
    strtitinf='$s$=Root mean squared error';
    annotation(h2,'textbox',dim,'FitBoxToText','on','String',strtitinf,'Interpreter','latex','FontSize',fs1);
    dimstrs = [posp .83 0.1 0.1];
    annotation(h2,'textbox',dimstrs,'FitBoxToText','on','String',strs,'Interpreter','latex','FontSize',fs);

    % Standard error section
    rr=0.03;
    dimstrs = [posp .70+rr 0.1 0.1];

    if intercept == true
        dimstderra = [posp .63+rr 0.1 0.1];
        annotation(h2,'textbox',dimstderra,'FitBoxToText','on','String',strstderra,'Interpreter','latex','FontSize',fs);
        strSE='Standard error of intercept and slope';
    else
        strSE='Standard error of slope';
    end
    dimstderrb = [posp .54+rr 0.1 0.1];
    annotation(h2,'textbox',dimstderrb,'FitBoxToText','on','String',strstderrb,'Interpreter','latex','FontSize',fs);

    annotation(h2,'textbox',dimstrs,'FitBoxToText','on','String',strSE,'Interpreter','latex','FontSize',fs1);


    % t statistic section
    dimstrtstat = [posp .41+rr 0.1 0.1];
    strtstat='$t$-statistic for a two-sided test with $H_0$: coefficient is zero';
    annotation(h2,'textbox',dimstrtstat,'FitBoxToText','on','String',strtstat,'Interpreter','latex','FontSize',fs1);

    if intercept == true
        strta=['\boldmath{$t_a= \frac{a}{s_a} = \frac{' num2str(a) '}{' num2str(stderra) '}= $ $' num2str(ta) '$}' '$\qquad p$-value=$Pr(|T({' num2str(lenx-2) '})|>|t_a|)=' num2str(pvala) '$'];
        dimstrta = [posp .34+rr 0.1 0.1];
        annotation(h2,'textbox',dimstrta,'FitBoxToText','on','String',strta,'Interpreter','latex','FontSize',fs);
    end

    strtb=['\boldmath{$t_b= \frac{b}{s_b} = \frac{' num2str(b) '}{' num2str(stderrb) '}= $ $' num2str(tb) '$}' '$\qquad p$-value=$Pr(|T({' num2str(lenx-2) ')}|>|t_b|)=' num2str(pvalb) '$'];
    dimstrtb = [posp .27+rr 0.1 0.1];
    annotation(h2,'textbox',dimstrtb,'FitBoxToText','on','String',strtb,'Interpreter','latex','FontSize',fs);

    % Confidence interval section
    dimstrconfint = [0.01 .28 0.1 0.1];
    strconfint='Confidence intervals';
    annotation(h2,'textbox',dimstrconfint,'FitBoxToText','on','String',strconfint,'Interpreter','latex','FontSize',fs1);

    % confint for a
    if intercept == true
        numstrafmt=num2str(a,'%.4f');
        numstrstderra=num2str(stderra,'%.4f');

        strconfinta=['\boldmath{$Pr(a+t_{' num2str(alpha/2) '}('   num2str(lenx-2)  ') \times s_a < \alpha <a+t_{' num2str(1-alpha/2) '}('  num2str(lenx-2)  ')\times s_a $ $)=' '$}' ...
            '\boldmath{$Pr(' numstrafmt num2str(tinvinf) '\times' numstrstderra  '<\alpha< '  numstrafmt '+' num2str(tinvsup) '\times' numstrstderra '$)= $' num2str(1-alpha) '$}'];
        dimstrconfinta = [0.005 .20 0.1 0.1];
        annotation(h2,'textbox',dimstrconfinta,'FitBoxToText','on','String',strconfinta,'Interpreter','latex','FontSize',fs);

        dimstrconfintafin = [0.30 .15 0.1 0.1];
        strconfintafin=['\boldmath{$Pr(' num2str(tinfa)  '< \alpha < '  num2str(tsupa) '$)= $' num2str(1-alpha) '$}'];
        annotation(h2,'textbox',dimstrconfintafin,'FitBoxToText','on','String',strconfintafin,'Interpreter','latex','FontSize',fs1);
    end

    % confint for b
    numstrbfmt=num2str(b,'%.4f');
    numstrstderrb=num2str(stderrb,'%.4f');
    strconfintb=['\boldmath{$Pr(b+t_{' num2str(alpha/2) '}('   num2str(lenx-2)  ') \times s_b < \beta <b+t_{' num2str(1-alpha/2) '}('  num2str(lenx-2)  ')\times s_b $ $)=' '$}' ...
        '\boldmath{$Pr(' numstrbfmt num2str(tinvinf) '\times' numstrstderrb  '<\beta< '  numstrbfmt '+' num2str(tinvsup) '\times' numstrstderrb '$)= $' num2str(1-alpha) '$}'];
    dimstrconfintb = [0.005 .05 0.1 0.1];
    annotation(h2,'textbox',dimstrconfintb,'FitBoxToText','on','String',strconfintb,'Interpreter','latex','FontSize',fs);

    dimstrconfintbfin = [0.30 .0001 0.1 0.1];
    strconfintbfin=['\boldmath{$Pr(' num2str(tinfb)  '< \beta < '  num2str(tsupb) '$)= $' num2str(1-alpha) '$}'];
    annotation(h2,'textbox',dimstrconfintbfin,'FitBoxToText','on','String',strconfintbfin,'Interpreter','latex','FontSize',fs1);

end

% Show plot of (x,y) coordinates, fitted line and residuals
if plots==true
    r2str=['(R2=' num2str(R2) ')'];
    if timeseries==true
        labx='t';
    else
        labx='x_i';
    end
    figure
    hold('on')

    if interpType==1
        plot(x,y,'o','LineWidth',2,'DisplayName','xyori')
        h=refline(b,a);
        %h=lsline;
        h.Color='r';

        if b>0
            abx=['(' num2str(a) '+' num2str(b) labx ')'];
        else
            abx=['(' num2str(a) num2str(b) labx ')'];
        end
        leg={['(' labx(1) ',y) coordinates'],['Regression line ' abx] ['Residuals ' r2str]};
        % plot([x x]',[yhat y]','k','DisplayName','Residuals')
        plot([x x]',[yhat y]','k','DisplayName',leg{3})
    end

    if interpType==2
        if b>0
            abx=['exp(' num2str(a) '+' num2str(b) labx ')'];
        else
            abx=['exp(' num2str(a) num2str(b) labx ')'];
        end

        leg={'(t,y_t) coordinates',['Exponential fit ' abx] ['Residuals ' r2str]};
        plot(x,yori,'o','LineWidth',2)
        plot(x,yhat,'DisplayName','fit')
        plot([x x]',[yhat yori]','k','DisplayName',leg{3})



    elseif interpType==3
        abx=[num2str(exp(a)) 't^{' num2str(b) '}'];
        leg={'(t,y_t) coordinates',['Power fit ' abx] ['Residuals ' r2str]};

        plot(xori,yori,'o','LineWidth',2)
        plot(xori,yhat,'r','DisplayName','fit')
        plot([xori xori]',[yhat yori]','k','DisplayName',leg{3})

    else
    end
    drawnow
    clickableMultiLegend(leg,'Location','best')
end

end
%FScategory:GUI