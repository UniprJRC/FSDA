function outFORE=forecastH(y,X,Z,varargin)
%forecastH produce forecasts with confidence bands for regression model with heteroskedasticity
%
%<a href="matlab: docsearchFS('forecastH')">Link to the help function</a>
%
%  Required input arguments:
%
%    y:         Response variable. Vector. Response variable, specified as
%               a vector of length n, where n is the number of
%               observations. Each entry in y is the response for the
%               corresponding row of X.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%  X :          Predictor variables in the regression equation. Matrix. Matrix of
%               explanatory variables (also called 'regressors') of
%               dimension n x (p-1) where p denotes the number of
%               explanatory variables including the intercept. Rows of X
%               represent observations, and columns represent variables. By
%               default, there is a constant term in the model, unless you
%               explicitly remove it using input option intercept, so do
%               not include a column of 1s in X. Missing values (NaN's) and
%               infinite values (Inf's) are allowed, since observations
%               (rows) with missing or infinite values will automatically
%               be excluded from the computations.
%     Z :       Predictor variables in the skedastic equation. Matrix. n x
%               r matrix or vector of length r. If Z is a n x r matrix it
%               contains the r variables which form the scedastic function.
%               If Z is a vector of length r it contains the indexes of
%               the columns of matrix X which form the scedastic function.
%               Therefore, if for example the explanatory variables
%               responsible for heteroscedisticity are columns 3 and 5
%               of matrix X, it is possible to use both the sintax
%                    regressH(y,X,X(:,[3 5]))
%               or the sintax
%                    regressH(y,X,[3 5])
%
%  Optional input arguments:
%
%   bsb :       list of units forming the initial
%               subset. Vector or []. If bsb=[] (default) then all units
%               are used to produce parameter estimates, else, just the
%               units forming bsb are used
%               Example - 'bsb',20:50
%               Data Types - double
%
%      conflev : confidence level for the confidence bands. Scalar.
%                A number between 0 and 1 which defines the confidence
%                level which is used to produce the bands. The default
%                value of conflev is 0.99.
%               Example - 'conflev',0.999
%               Data Types - double
%
% intercept :   Indicator for constant term. true (default) | false.
%               Indicator for the constant term (intercept) in the fit,
%               specified as the comma-separated pair consisting of
%               'Intercept' and either true to include or false to remove
%               the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
%    outH:      output from fitted heteroskedatic model. Struct.
%               it is possible to supply the output produced by
%               functions, regressH or regressHart or regressHhar or FSRHeda.
%               Note that if input optional argument outH is supplied the model is not fitted
%               and the parameter estimates are taken from outH.
%                 Example - out=regress(y,X,Z); 'outH',out
%                 Data Types - struct
%
% originalScale: confidence band in original or transformed scale. Boolean.
%               If originalScale is true, plot is shown in the original
%               scale (default). If originalScale is false forecasts are
%               shown on the transformed scale.
%                 Example - 'originalScale',false
%                 Data Types - boolean
%
%    selcolX:   column of matrix X for which confidence band is wanted.
%               Integer in the set $1, 2, \ldots, p-1$.
%               Scalar which identifies the column of X
%               to put in x axis of the plot. Default value of selcolX is
%               1.
%                 Example - 'selcolX',2
%                 Data Types - double
%
%   typeH:      Parametric function to be used in the skedastic equation.
%               Character or string.
%               If typeH is 'art' (default) than the skedastic function is
%               modelled as follows
%               \[
%               \sigma^2_i = \sigma^2 (1 + \exp(\gamma_0 + \gamma_1 Z(i,1) +
%                           \cdots + \gamma_{r} Z(i,r)))
%               \]
%               on the other hand, if typeH is 'har' then traditional
%               formulation due to Harvey is used as follows
%               \[
%               \sigma^2_i = \exp(\gamma_0 + \gamma_1 Z(i,1) + \cdots +
%                           \gamma_{r} Z(i,r)) =\sigma^2 (\exp(\gamma_1
%                           Z(i,1) + \cdots + \gamma_{r} Z(i,r))
%               \]
%               Remark. Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%               Example - 'typeH','har'
%               Data Types - character or string
%
%
%  Output:
%
%         outFORE:   structure which contains the following fields
%                    outFORE.conf = matrix of size nx4 referred to
%                                original space.
%                                1st column = ordered X values;
%                                2nd column = fitted values;
%                                3rd column = lower confidence band;
%                                4th column = upper confidence band.
%                   outFORE.confW = matrix of size nx4 referred to
%                                transformed space.
%                                1st column = ordered X values;
%                                2nd column = fitted values;
%                                3rd column = lower confidence band;
%                                4th column = upper confidence band.
%
%
% See also regressH, regressHart, regressHhar, FSRHeda
%
% References:
%
% Greene, W.H. (1987), "Econometric Analysis", Prentice Hall. [5th edition,
% section 11.7.1 pp. 232-235, 7th edition, section  9.7.1 pp. 280-282]
%
% Atkinson, A.C., Riani, M. and Torti, F. (2016), Robust methods for
% heteroskedastic regression, "Computational Statistics and Data Analysis",
% Vol. 104, pp. 209-222, http://dx.doi.org/10.1016/j.csda.2016.07.002 [ART]
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('forecastH')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % forecastH with all default options.
    close all
    load tradeH.mat
    y=tradeH{:,2};
    X=tradeH{:,1};
    X=X./max(X);
    Z=log(X);
    fore=forecastH(y,X,Z)
%}

%{
    %%  Example of use of option bsb.
    close all
    load tradeH.mat
    y=tradeH{:,2};
    X=tradeH{:,1};
    X=X./max(X);
    Z=log(X);
    n=length(y);
    bsb=1:n;
    % outl = the outliers
    outl=[225   660];
    bsb(outl)=[];
    % call of forecastH with option bsb
    fore=forecastH(y,X,Z,'bsb',bsb);
%}

%{
    %%  Example of use of option typeH.
    close all
    load tradeH.mat
    y=tradeH{:,2};
    X=tradeH{:,1};
    X=X./max(X);
    Z=log(X);
    % Use Harvery's parametrization
    fore=forecastH(y,X,Z,'typeH','har');
%}

%{
    % Example of use of option outH.
    close all
    load tradeH.mat
    y=tradeH{:,2};
    X=tradeH{:,1};
    X=X./max(X);
    Z=log(X);
    outEDA=FSRHeda(y,X,Z,0,'init',round(length(y)/2));
    % use the output of outEDA
    fore=forecastH(y,X,Z,'outH',outEDA);
%}

%{
    %% Monthly credit card expenditure for 100 individuals.
    % Example of use of option selcolX
    load('TableF91_Greene');
    data=TableF91_Greene{:,:};
    n=size(data,1);
    
    % Linear regression of monthly expenditure on a constant, age, income
    % its square and a dummy variable for home ownership using the 72
    % observations for which expenditure was nonzero produces the residuals
    % plotted below
    
    X=zeros(n,4);
    X(:,1)=data(:,3);%age
    X(:,2)=data(:,6);% Own rent (dummy variable)
    X(:,3)=data(:,4);% Income
    X(:,4)=(data(:,4)).^2; %Income  square
    y=data(:,5); % Monthly expenditure
    
    % Select the 72 observations for which expenditure was nonzero
    sel=y>0;
    X=X(sel,:);
    y=y(sel);
    close all
    disp('Multiplicative Heteroskedasticity Model')
    % Plot of forecasts against column 4
    warning('off')
    out=forecastH(y,X,[3 4],'typeH','har','selcolX',4);     
    warning('on')
%}

%% Beginning of code

bsb=[];
conflev=0.95;
n=length(y);
typeH='art';
outH='';
intercept=true;
originalScale=true;
selcolX=1;
if size(Z,1)~=n
    Z=[ones(n,1) X(:,Z)];
end

if nargin>3
    options=struct('intercept',intercept, ...
        'outH',outH,'bsb',bsb,'conflev',conflev,'typeH',typeH, ...
        'originalScale',true,'selcolX',selcolX);

    UserOptions=varargin(1:2:length(varargin));
    % Check if user options are valid options
    chkoptions(options,UserOptions)


    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end
    bsb=options.bsb;
    typeH=options.typeH;
    conflev=options.conflev;
    intercept=options.intercept;
    % Show plot on original or transformed scale
    originalScale=options.originalScale;
    % Column of X to show on the x axis of the plot
    selcolX=options.selcolX;

    % check if input option outH exists
    chkoutH = 2*find(strcmpi('outH',UserOptions), 1);
    if ~isempty(chkoutH)
        % In this case the heteroskedastic model has already been fitted
        out=options.outH;
        bsb=options.bsb;
        % Get the type of heteroskedasticity which has been used
        typeH=out.typeH;
    else

        chkbsb = 2*find(strcmpi('bsb',UserOptions));
        chktypeH = 2*find(strcmpi('typeH',UserOptions));
        chkconflev=2*find(strcmpi('conflev',UserOptions));
        chkoriscale=2*find(strcmpi('originalScale',UserOptions));
        chkselcolX=2*find(strcmpi('selcolX',UserOptions));

        % remove from varargin1 options typeH, bsb, conflev,
        % selcolX and originalScale if they exist
        varargin1=varargin;
        varargin1([chktypeH-1 chktypeH chkbsb-1 chkbsb ...
            chkconflev-1 chkconflev chkoriscale-1 chkoriscale ...
            chkselcolX-1 chkselcolX])=[];

        if isempty(options.bsb)
            yb=y;
            Xb=X;
            Zb=Z;
        else
            yb=y(bsb);
            Xb=X(bsb,:);
            Zb=Z(bsb,:);
        end

        if ~isempty(chktypeH) && strcmp(typeH,'har') ==1
            out=regressHhar(yb,Xb,Zb,varargin1{:});
        else
            out=regressHart(yb,Xb,Zb,varargin1{:});
        end
    end
else % In this case the heteroskedastic model has not been fitted yet
    % and bsb has not been specified
    out=regressHart(y,X,Z);
end

seq=1:n;

[~,indori]=sort(X(:,selcolX));

quant=norminv(1-(1-conflev)/2);


if isempty(bsb)
    bsb=seq;
end

outl=setdiff(seq,bsb);
Xout=X(outl,:);
yout=y(outl);

bsb=sort(bsb(:));

% Now produce the confidence bands
if isfield(out,'class') && strcmp(out.class,'FSRHeda') % Case 1 class is FSRHeda


    % Find column to extract from out.BB
    for j=1:size(out.BB,2)
        if isequal(bsb,rmmissing(out.BB(:,end-j+1)))
            break
        end
    end
    colToextract=size(out.BB,2)-j+1;

    sigma2=out.S2(colToextract,2);
    % inv of sqrt of weights (vector of length(n)
    sqweights = out.WEI(:,colToextract);
else % in this case the output comes from regressH, regressHart or regressHhar
    if strcmp(typeH,'art')
        omegahat = 1+exp(out.Gamma(1,1))*exp(Z*out.Gamma(2:end,1));
    else
        omegahat = exp(Z*out.Gamma(2:end,1));
    end
    sqweights = omegahat.^(-0.5);
    sigma2=out.sigma2;

end

bool=true(n,1);

if ~isempty(outl)
    bool(outl)=false;
    XoutW=Xout.*sqweights(outl);
    youtW=yout.*sqweights(outl);
end

Xgood=X(bsb,:);
ygood=y(bool);
sqweightsgood=sqweights(bool);
XgoodW=Xgood.*sqweightsgood;
ygoodW=ygood.*sqweightsgood;

if intercept==1
    Xfull=[ones(n,1) X];
else
    Xfull=X;
end

% XfullW and yW = X and y in the transformed scale
XW=X.*sqweights;
XfullW = Xfull.* sqweights;
yW = y .* sqweights;

% XXw=(X'*W*X)^-1 based on the units forming subset (without the outliers)
XXw=inv(XfullW(bool,:)'*XfullW(bool,:));

% Beta= Estimate of beta without the outliers
Beta=XfullW(bool,:)\yW(bool);


%  yhatW = fitted values in the transformed scale
yhatW=XfullW*Beta;
% yhat = fitted values in the original scale
yhat=Xfull*Beta;

hold('on')


% MSEt= sqrt of variance of prediction for each value X
% in the transformed scale
MSEt=zeros(n,1);
for i=1:n
    MSEt(i)=sqrt(sigma2.*(1+XfullW(i,:)*XXw*(XfullW(i,:)'))); %#ok<MINV>
end

% Fitted values in the original scale
Xord=X(indori,selcolX);
yhatord=yhat(indori);
% confidence bands in the original scale
lowconf=(yhatW-quant*MSEt)./sqweights;
lowconford=lowconf(indori);
upconf=(yhatW+quant*MSEt)./sqweights;
upconford=upconf(indori);

% Fitted values in the transformed scale
XWord=XW(indori,selcolX);
yhatordW=yhatW(indori);
% confidence bands in the transformed scale
lowconfW=(yhatW-quant*MSEt);
lowconfordW=lowconfW(indori);
upconfW=(yhatW+quant*MSEt);
upconfordW=upconfW(indori);


if originalScale==true
    % Plot the regression line in the original scale
    plot(Xord,yhatord)
    % plot good observations
    plot(Xgood(:,selcolX),ygood,'ro')
    % plot the outliers
    if ~isempty(outl)
        plot(Xout,yout,'kx','MarkerSize',12)
    end

    % Plot the upper and lower confidence bands (in the original scale)
    plot(Xord,lowconford,'k')
    plot(Xord,upconford,'k')

else
    plot(XWord,yhatordW)
    hold('on')
    plot(XgoodW,ygoodW,'ro','MarkerSize',8)
    if ~isempty(outl)
        plot(XoutW,youtW,'kx','MarkerSize',8)
    end
    % Plot the upper and lower confidence bands (in the transformed scale)
    plot(XWord,lowconfordW,'k')
    plot(XWord,upconfordW,'k')
end
xlabel(['Variable ' num2str(selcolX)])
hold('off')

% Store the fitted values and confidence bands

outFORE=struct;
outFORE.conf=[Xord yhatord lowconford upconford];
outFORE.confW=[XWord yhatordW lowconfordW upconfordW];
end
%FScategory:REG-Hetero