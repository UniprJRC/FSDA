function [outSC]=ScoreYJmle(y,X,varargin)
%ScoreYJmle computes the likelihood ratio test fof H_0=lambdaP=lambdaP0 and lambdaN=lambdaN0
%
%
%<a href="matlab: docsearchFS('ScoreYJmle')">Link to the help function</a>
%
%
% The transformations for negative and positive responses were determined
% by Yeo and Johnson (2000) by imposing the smoothness condition that the
% second derivative of zYJ(λ) with respect to y be smooth at y = 0. However
% some authors, for example Weisberg (2005), query the physical
% interpretability of this constraint which is oftern violated in data
% analysis. Accordingly, Atkinson et al (2019) and (2020) extend the
% Yeo-Johnson transformation to allow two values of the transformations
% parameter: λN for negative observations and λP for non-negative ones.
% $\lambda$ is the transformation parameter (scalar) for all the
% obseravtions (positive adn negative).
% $\lambda_P$ is the transformation parameter for positive observations.
% $\lambda_N$ is the transformation parameter for negative observations.
% SSR is the residual sum of squares of the model which regresses
% $z(λ)$ against X.
% SSF is the residual sum of squares of the model which regresses
% $z(\hat λ_{MLE})$ against $X$ where $\lambda_{MLE}$ is the vector
% of length 2 which contains the MLE of $\lambda_P$ and $\lambda_N$
% ScoreYJmle computes Num/Den where Num and Den are defined as follows:
% Num=(SSR-SSF)/2  and Den=SSF/(n-p-2) where p is the number of columns of
% matrix X (including intercept).
%
%
%
%
%  Required input arguments:
%
%    y:         Response variable. Vector. A vector with n elements that
%               contains the response
%               variable.  It can be either a row or a column vector.
%    X :        Predictor variables. Matrix. Data matrix of explanatory
%               variables (also called 'regressors')
%               of dimension (n x p-1). Rows of X represent observations, and
%               columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with missing or infinite values will
%               automatically be excluded from the computations.
%
%  Optional input arguments:
%
%    intercept :  Indicator for constant term. true (default) | false.
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as the comma-separated pair consisting of
%                 'Intercept' and either true to include or false to remove
%                 the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
%        la  :  transformation parameter. Vector. It specifies for which
%               values of the transformation parameter it is necessary to
%               compute the score test. Default value of lambda is la=[-1
%               -0.5 0 0.5 1]; that is the five most common values of
%               lambda
%               Example - 'la',[0 0.5]
%               Data Types - double
%
%    usefmin :  use solver to find MLE of lambda. Boolean or struct.
%               if usefmin is true or usefmin is a struct it is
%               possible to use MATLAB solvers fminsearch or fminunc to
%               find the maximum likelihood estimates of $\lambda_P$ and
%               $\lambda_N$. The default value of usefmin is false that is
%               solver is not used and the likelihood is evaluated at the
%               grid of points with steps 0.01.
%               If usefmin is a structure it may contain the following
%               fields:
%               usefmin.MaxIter = Maximum number of iterations (default is 1000).
%               usefmin.TolX   = Termination tolerance for the parameters
%                   (default is 1e-7).
%               usefmin.solver = name of the solver. Possible values are
%                   'fminsearch' (default) and 'fminunc'. fminunc needs the
%                   optimization toolbox.
%               usefmin.displayLevel = amount of information displayed by
%                   the algorithm. possible values are 'off' (displays no
%                   information, this is the default), 'final' (displays
%                   just the final output) and 'iter' (displays iterative
%                   output to the command window).
%               Example - 'usefmin',true
%               Data Types - boolean or struct
%
% sseReducedModel: sum of squares of residuals of reduced model. Vector.
%               Vector with the same length of input vector lambda
%               containing the sum of squares of residuals of the reduced
%               model. The default value of sseReducedModel is an empty
%               value that is this quantity is computed by this routine.
%               Example - 'sseReducedModel',[20.2 30.3 12.8]
%               Data Types - empty value or double
%
%       nocheck : Check input arguments. Scalar.
%               If nocheck is equal to 1 no check is performed on
%                 matrix y and matrix X. Notice that y and X are left
%                 unchanged. In other words the additional column of ones
%                 for the intercept is not added. As default nocheck=0.
%               Example - 'nocheck',1
%               Data Types - double
%
%  Output:
%
%  The output consists of a structure 'outSC' containing the following fields:
%
%        outSC.Score =       score tests. Vector.
%                            Column vector of length(la) which
%                            contains the value of the likelihood ratio test for each
%                            value of lambda specified in optional input
%                            parameter la.
%        outSC.laMLE =       score tests. Matrix.
%                            Matrix of dimension length(la)-by-2 which
%                            contains the value of maximum likelihood
%                            estimate of $\lambda_P$ and $\lambda_N$
%                            for each value of lambda specified in optional
%                            input parameter la. First column refers to
%                            $\lambda_P$ and second column to $\lambda_N$
%
% See also: FSRfan, Score, ScoreYJ, ScoreYJpn,
%
% References:
%
% Yeo, I.K. and Johnson, R. (2000), A new family of power
% transformations to improve normality or symmetry, "Biometrika", Vol. 87,
% pp. 954-959.
% Atkinson, A.C. Riani, M., Corbellini A. (2019), The analysis of
% transformations for profit-and-loss data, Journal of the Royal
% Statistical Society, Series C, "Applied Statistics",
% https://doi.org/10.1111/rssc.12389
% Atkinson, A.C. Riani, M. and Corbellini A. (2020), The Box-Cox
% Transformation: Review and Extensions, "Statistical Science", in press.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('ScoreYJmle')">Link to the help function</a>
%
%$LastChangedDate:: 2017-11-17 15:01:40 #$: Date of the last commit

% Examples


%{
    rng('default')
    rng(1)
    n=100;
    yori=randn(n,1);
    % Transform the value to find out if we can recover the true value of
    % the transformation parameter
    la=0.5;
    ytra=normYJ(yori,[],la,'inverse',true);
    % Start the analysis
    X=ones(n,1);
    [outSCmle]=ScoreYJmle(ytra,X,'intercept',0);

    la=[-1 -0.5 0 0.5 1]';
    Comb=[la outSCmle.Score(:,1)];
    CombT=array2table(Comb,'VariableNames',{'la', 'FtestLR'});
    disp(CombT)
%}

%{
    %% Ex in which positive and negative observations require different lambdas.
    rng(1000)
    n=100;
    y=randn(n,1);
    % Transform in a different way positive and negative values
    lapos=0;
    ytrapos=normYJ(y(y>=0),[],lapos,'inverse',true);
    laneg=1;
    ytraneg=normYJ(y(y<0),[],laneg,'inverse',true);
    ytra=[ytrapos; ytraneg];

    % Start the analysis
    X=ones(n,1);
    la=[-1:0.25:1]';
    [outSCmle]=ScoreYJmle(ytra,X,'intercept',0,'la',la);
    Pval=fcdf(outSCmle.Score,2,n-2,'upper');
    Comb=[la outSCmle.Score(:,1) Pval];
    
    CombT=array2table(Comb,'VariableNames',{'la','FtestLR' 'Pvalues'});
    disp(CombT)
    disp('The test is significant for all values of lambda')
    disp('This may indicate that the data need two separate lambdas for pos and neg observations')
%}



%% Beginning of code


nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

la=[-1 -0.5 0 0.5 1];

usefmin=true;
sseReducedModel=[];

if nargin>2
    
    options=struct('la',la,'nocheck',0,'intercept',0,...
        'usefmin',usefmin,'sseReducedModel',sseReducedModel);
    
    UserOptions=varargin(1:2:length(varargin));
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:ScoreYJmle:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    
    % Write in structure 'options' the options chosen by the user
    if nargin > 2
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
    end
    
    sseReducedModel=options.sseReducedModel;
    la=options.la;
    usefmin=options.usefmin;
end


%% Find MLE of la_P and la_N once and for all

if isstruct(usefmin) || (islogical(usefmin) && usefmin==true)
    LA=boxcoxR(y,X,'family','YJpn','plots',0,'nocheck',1,...
        'usefmin',usefmin);
    
else
    lamax=2;
    step=0.02;
    LA=boxcoxR(y,X,'family','YJpn','plots',0,'nocheck',1,...
        'laseqPos',-lamax:step:lamax,'laseqNeg',-lamax:step:lamax);
end

% Transform the data using MLE
zMLE=normYJpn(y,[],LA.lahat,'Jacobian',true);

% Find sum of squares of residuals using zMLE
beta=X\zMLE;
residualsF = zMLE - X*beta;

% Sum of squares of residuals for FULL model
SSeF = norm(residualsF)^2;

%% Loop over the values of lambda_0

%  LikRatio= vector of length  lla which contains the likelihood ratio test
%  for the values of \lambda specified in vector la
lla=length(la);
LikRatio=NaN(lla,1);

% Initialize things which are constant across the different values of lambda
nonnegs = y >= 0;
negs = ~nonnegs;
ynonnegs=y(nonnegs);
ynegs=y(negs);

logG=sum(  sign(y) .* log(abs(y)+1)   )/n;
vneg=-ynegs+1;
vpos=ynonnegs+1;
logvpos=log(vpos);
logvneg=log(vneg);
G=exp(logG);

% loop over the values of \lambda
for i=1:lla
    
    if isempty(sseReducedModel)
        z=y; % Initialized z and w
        lai=la(i);
        Glaminus1=G^(lai-1);
        q=lai*Glaminus1;
        twomlambdai=2-lai;
        
        % Compute transformed values under the null hypotesis
        % transformation for non negative values
        if abs(lai)>1e-8  % if la is different from 0
            % vposlai=vpos.^lai;
            vposlai=exp(lai*logvpos);
            znonnegs=(vposlai-1)/q;
            z(nonnegs)=znonnegs;
        else % if la is equal to 0
            znonnegs=G*logvpos;
            z(nonnegs)=znonnegs;
        end
        
        % Transformation and constructed variables for negative values
        if   abs(lai-2)>1e-8 % la not equal 2
            % vnegtwomlambdai=vneg.^twomlambdai;
            vnegtwomlambdai=exp(twomlambdai*logvneg);
            qneg=twomlambdai* Glaminus1;
            znegs=(1-vnegtwomlambdai )  /qneg;
            z(negs)=znegs;
        else  % la equals 2
            znegs=-logvneg/G;
            z(negs)=znegs;
        end
        
        % Compute residual sum of squares for null (reduced) model
        betaR=X\z;
        residualsR = z - X*betaR;
        % Sum of squares of residuals
        SSeR = norm(residualsR)^2;
    else
        SSeR=sseReducedModel(i);
    end
    
    Ftestnum=(SSeR-SSeF)/2;
    Ftestden=SSeF/(n-p-2);
    Ftest=Ftestnum/Ftestden;
    LikRatio(i,1)=Ftest;
    
end

% Store values of the score test inside structure outSC
outSC.Score=LikRatio;
% Store values of the MLE of lambda inside structure outSC
outSC.laMLE=LA.lahat;
end
%FScategory:REG-Transformations