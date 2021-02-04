function [out]=addt(y,X,w,varargin)
%addt produces the t test for an additional explanatory variable
%
%<a href="matlab: docsearchFS('addt')">Link to the help page for this function</a>
%
% Required input arguments:
%
%       y:  A vector with n elements that contains the response variable.
%           y can be both a row of column vector.
%       X:  Data matrix of explanatory variables (also called
%           'regressors').
%           Rows of X represent observations and columns represent
%           variables.
%           Missing values (NaN's) and infinite values (Inf's) are allowed,
%           since observations (rows) with missing or infinite values will
%           automatically be excluded from the computations.
%       w:  added variable. Vector. n-x-1 vector containing the additional
%           explanatory variable whose t test must be computed.
%
% Optional input arguments:
%
%    intercept :  Indicator for constant term. true (default) | false.
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as the comma-separated pair consisting of
%                 'Intercept' and either true to include or false to remove
%                 the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
%   la:         Transformation parameter. Scalar | '' (empty value).
%               It specifies for which Box Cox
%               transformation parameter it is necessary to compute the t
%               statistic for the additional variable. If la is an empty
%               value (default) no transformation is used.
%               Example - 'la',0.5 tests square root transformation
%               Data Types - double
%
%   plots:      Plot on the screen. Scalar.
%               If plots=1 the added variable
%               plot is produced else (default) no plot is produced.
%               Example - 'plots',1
%               Data Types - double
%
%   units:      Units to remove. Vector.
%               Vector containing the list of
%               units which has to be removed in the computation of the
%               test. The default is to use all units
%               Example - 'units',[1,3] removes units 1 and 3
%               Data Types - double
%
%   textlab:    Labels of units in the plot. Boolean. If textlab=false
%               (default) no text label is written on the plot
%               for units else text label of units are added on the plot
%               Example - 'textlab',true
%               Data Types - boolean
%
%   FontSize:   Label font size inside plot. Scalar. It controls the
%               fontsize of the labels of the axes and eventual plot
%               labels. Default value is 10
%               Example - 'FontSize',14
%               Data Types - double
%
%   SizeAxesNum: Font size of axes numbers. Scalar. It controls the
%               fontsize of the numbers of the
%               axes. Default value is 10
%               Example - SizeAxesNum,12
%               Data Types - double
%
% nocheck :       Check input arguments. Scalar.
%               If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additional column of ones
%               for the intercept is not added. As default nocheck=0.
%               Example - 'nocheck',1
%               Data Types - double
%
% Output:
%
%         out:   structure which contains the following fields
%
%       out.b=          estimate of the slope for additional explanatory
%                       variable
%       out.S2add=  estimate of $s^2$ of the model which contains the
%                       additional explanatory variable
%       out.Tadd=         t statistic for additional explanatory variable
%       out.pval=         p-value of the t statistic
%
%
% See also FSRaddt
%
% References:
%
% Atkinson, A.C. and Riani, M. (2000), "Robust Diagnostic Regression
% Analysis", Springer Verlag, New York.
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('addt')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% addt with all default options.
    % Compute the t test for an additional regressor.
    XX=load('wool.txt');
    y=log(XX(:,end));
    X=XX(:,1:end-2);
    w=XX(:,end-1);
    [out]=addt(y,X,w);
    
    % out.Tadd (equal to -8.9707) is exactly equal to stats.tstat.t(4)
    % obtained as
    
    whichstats = {'tstat','mse'};
    stats = regstats(y,XX(:,1:end-1),'linear',whichstats);
    
    % Similarly out.S2add (equal to 0.0345) is exactly equal to stats.mse (estimate of
    % \sigma^2 for augmented model)
%}

%{
    %% addt with optional arguments.
    % Excluding one observation from the sample; compare the added variable plot
    % based on all units with that which excludes unit 43.
    load('multiple_regression.txt');
    y=multiple_regression(:,4);
    X=multiple_regression(:,1:3);
    [out]=addt(y,X(:,2:3),X(:,1),'plots',1,'units',[43],'textlab',true);
%}

%{
    %% Excluding more than one observation from the sample.
    % Compare the added variable plot based on all units with that which excludes units
    % 9,21,30,31,38 and 47.
    load('multiple_regression.txt');
    y=multiple_regression(:,4);
    X=multiple_regression(:,1:3);
    [out]=addt(y,X(:,2:3),X(:,1),'plots',1,'units',[9 21 30 31 38 47]','textlab',true);
%}
%
%
%% Beginning of code

% User options

if nargin<3
    error('FSDA:addt:missingInputs','A required input argument is missing.')
end

la=1;
units=[]';
textlab=false;
FontSize=10;
SizeAxesNum=10;
nocheck=0;
plots=0;


if nargin > 3
    
    if coder.target('MATLAB')
        % the 'options' struct in this implementation is provided by the for loop
        options=struct('intercept',1,'plots',plots,'la',la,'units',units,'textlab',textlab,...
            'FontSize',FontSize,'SizeAxesNum',SizeAxesNum,'nocheck',nocheck);
    end
    
    for i=1:2:(length(varargin)-1)
        options.(varargin{i})=varargin{i+1};
    end
    
    
    la=options.la;
    plots=options.plots;
    units=options.units;
    textlab=options.textlab;
    
    % FontSize = font size of the axes labels
    FontSize =options.FontSize;
    % FontSizeAxes = font size for the axes numbers
    SizeAxesNum=options.SizeAxesNum;
end
%% t test for an additional explanatory variable

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

[~, R] = qr(X,0);
E = X/R;
A = -E*E';
sel=1:n;
siz = size(A);
% Find linear indexes
% It is better to compute linind directly rather than calling sub2ind
% linind=sub2ind(siz,sel,sel);
linind = sel + (sel - 1).*siz(1);


A(linind)=1+A(linind);
% Notice that:
% -E*E' = matrix -H = -X*inv(X'X)*X' computed through qr decomposition
% A = Matrix I - H


if ~isempty(la)
    la1=la(1);
    %geometric mean of the y
    G=exp(mean(log(y)))+ 0i; %G is complex;
    %  if la1==0
    if la==0
        z=G*log(y);
    else
        z=(y.^la1-1)/(la*G^(la1-1));
    end
else
    z=y;
end

Az=A*z;
r=z'*Az;
Aw=A*w;
zAw=z'*Aw;
wAw=w'*Aw;

if wAw <1e-12
    Sz_square=NaN;
    Tl=NaN;
    b=NaN;
    pval=NaN;
    if coder.target('MATLAB')
        warning('FSDA:addt:NearlySingularMatrix','The augmented X matrix is nearly singular');
    end
else
    % b=regress(Az,Aw);
    b=zAw/wAw;
    
    Sz=sqrt(r-zAw^2/wAw); % See Atkinson (1985) p. 98
    Sz_square=Sz^2/(n-p-1);
    
    if abs(real(Sz)) > 0.0000001
        % Compute t-statistic
        Tl=zAw*sqrt(n-p-1)/(Sz*sqrt(wAw));
        % Compute p-value of t-statistic
        pval=2*(1-tcdf(abs(Tl),n-p-1));
    else
        Tl=NaN;
        pval=NaN;
    end
end


% Store results in structure out.
out.b=b;
out.S2add=Sz_square;
out.Tadd=Tl;
out.pval=pval;


%% Added variable plot

if coder.target('MATLAB') && plots==1
    if ~isempty(units)
        sel=setdiff(1:length(y),units);
        [outsel]=addt(y(sel),X(sel,2:end),w(sel),'plots',0);
        
        plot(Aw(sel),Az(sel),'+b','MarkerSize',FontSize);
        hold('on');
        plot(Aw(units),Az(units),'or','MarkerSize',6,'MarkerFaceColor','r');
        
        xlimits = get(gca,'Xlim');
        % Superimpose line based on all units
        line(xlimits , b.*xlimits,'Color','r','LineWidth',2);
        % Superimposed line based on reduced set of units
        line(xlimits , outsel.b.*xlimits,'LineWidth',2);
        if textlab==true
            text(Aw(units)+0.05,Az(units),num2str(units),'FontSize',FontSize);
        end
        set(gca,'FontSize',SizeAxesNum)
        
        xlabel('Aw','Fontsize',FontSize);
        ylabel('Ay','Fontsize',FontSize);
        % Format for the legend
        forleg='%11.3g';
        forleg1='%3.2g';
        
        legend('Normal units','Excluded units',['Fit on all units tstat=' num2str(Tl,forleg) ' (pval=' num2str(pval,forleg1) ')'],...
            ['Fit on subset tstat=' num2str(outsel.Tadd,forleg) ' (pval=' num2str(outsel.pval,forleg1) ')'])
        hold('off');
        
        %olsline(2)
    else
        plot(Aw,Az,'+');
        xlimits = get(gca,'Xlim');
        % Superimpose line based on all units
        line(xlimits, b.*xlimits,'LineWidth',3);
        set(gca,'FontSize',SizeAxesNum)
        
        xlabel('Aw','Fontsize',FontSize);
        ylabel('Ay','Fontsize',FontSize);
        
    end
    
    % Make the legends clickable.
    hLines = findobj(gca, 'type', 'line');
    eLegend = cell(length(hLines), 1);
    for iLines = 1:length(hLines)
        eLegend{iLines} = get(hLines(iLines), 'DisplayName');
    end
    clickableMultiLegend(hLines, eLegend{:});
    
    % and freeze the scaling at the current limits
    axis(axis); % equivalent to "axis manual";
    
end

end
%FScategory:REG-Regression