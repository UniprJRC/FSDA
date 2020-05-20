function [outSC]=ScoreYJall(y,X,varargin)
%Computes all the 4 score tests for YJ transformation
%
%<a href="matlab: docsearchFS('ScoreYJall')">Link to the help function</a>
%
% The transformations for negative and positive responses were determined
% by Yeo and Johnson (2000) by imposing the smoothness condition that the
% second derivative of zYJ(λ) with respect to y be smooth at y = 0. However
% some authors, for example Weisberg (2005), query the physical
% interpretability of this constraint which is oftern violated in data
% analysis. Accordingly, Atkinson et al (2019) and (2020) extend the
% Yeo-Johnson transformation to allow two values of the transformations
% parameter: λN for negative observations and λP for non-negative ones.
% ScoreYJall computes:
% 1) the global t test associated with the constructed variable for
% λ=λP=λN.
% 2) the t test associated with the constructed variable computed assuming
% a different transformation for positive observations keeping the value of
% the transformation parameter for negative observations fixed. In short we
% call this test, "test for positive observations".
% 3) the t test associated with the constructed variable computed assuming
% a different transformation for negative observations keeping the value of
% the transformation parameter for positive observations fixed. In short we
% call this test, "test for negative observations".
% 4) the F test for the joint presence of the two constructed variables
% described in points 2) and 3.
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
%      scoremle: likelihood ratio test for the two different transformation
%                parameters $\lambda_P$ and $\lambda_N$. Boolean.
%                if scoremle is true it is possible to compute the
%                likelihood ratio test. In this case the residual sum of
%                squares of the null model bsaed on a single trasnformation
%                parameter $\lambda$ is compared with the residual sum of
%                squares of the model based on data transformed data using
%                MLE of $\lambda_P$ and $\lambda_N$. If scoremle is true it
%                is possible through following option usefmin, to control
%                the parameters of the optmization routine.
%               Example - 'scoremle',true
%               Data Types - logical
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
%        outSC.Score =       score tests. Matrix.
%                            Matrix of size length(la)-by-5 which
%                            contains the value of the score test for each
%                            value of lambda specified in optional input
%                            parameter la. The first column refers to the
%                            global test, the second to the test for
%                            positive observations, the third refers to the
%                            test for negative observations and the fourth
%                            column refers to the F test for the joint
%                            presence of the two constructed variables.
%                            If input option scoremle is true the fifth
%                            column will contain the exact likelihod ratio
%                            test based on the maximum likelihood estimates
%                            of the $\lambda_P$ and $\lambda_N$.
%                            If la is not specified, the number of rows of
%                            outSc.Score is equal to 5 and will contain the
%                            values of the score tests for the 5 most common
%                            values of lambda.
%        outSC.laMLE =       MLE of lambda. Vector.
%                            Vector of dimension 2 which
%                            contains the value of maximum likelihood
%                            estimate of $\lambda_P$ and $\lambda_N$.  This
%                            output is present only if input option
%                            scoremle is true.
%                            
%
% See also: FSRfan, Score, ScoreYJ, ScoreYjpn, fanBIC
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
%<a href="matlab: docsearchFS('ScoreYJall')">Link to the help function</a>
%
%$LastChangedDate:: 2017-11-17 15:01:40 #$: Date of the last commit

% Examples


%{
    %% Ex in which positive and negative observations require the same lambda.
    rng('default')
    rng(1)
    n=100;
    y=randn(n,1);
    % Transform the value to find out if we can recover the true value of
    % the transformation parameter
    la=0.5;
    ytra=normYJ(y,[],la,'inverse',true);
    % Start the analysis
    X=ones(n,1);
    [outSC]=ScoreYJ(ytra,X,'intercept',0);
    [outSCpn]=ScoreYJpn(ytra,X,'intercept',0);
    la=[-1 -0.5 0 0.5 1]';
    disp([la outSCpn.Score(:,1) outSC.Score outSCpn.Score(:,2)])
    % Comment: if we consider the 5 most common values of lambda
    % the value of the score test when lambda=0.5 is the only one which is not
    % significant. Both values of the score test for positive and negative
    % observations confirm that this value of the transformation parameter is
    % OK for both sides of the distribution.
%}

%{
    %% Ex in which positive and negative observation require different lambdas.
    rng(1000)
    n=100;
    y=randn(n,1);
    % Tranform in a different way positive and negative values
    lapos=0;
    ytrapos=normYJ(y(y>=0),[],lapos,'inverse',true);
    laneg=1;
    ytraneg=normYJ(y(y<0),[],laneg,'inverse',true);
    ytra=[ytrapos; ytraneg];

    % Start the analysis
    X=ones(n,1);
    [outSC]=ScoreYJ(ytra,X,'intercept',0);
    [outSCpn]=ScoreYJpn(ytra,X,'intercept',0);
    la=[-1 -0.5 0 0.5 1]';
    disp([la outSCpn.Score(:,1) outSC.Score outSCpn.Score(:,2)])
    % Comment: if we consider the 5 most common values of lambda
    % the value of the score test when lambda=0.5 is the only one which is not
    % significant. However when lambda=0.5 the score test for negative
    % observations is highly significant. The difference between the test for
    % positive and the test for negative is 2.7597+0.7744=3.5341, which is very
    % large. This indicates that the two tails need a different value of the
    % transformation parameter.
%}

%{
    % Extended score with all default options for the wool data.
    % Load the wool data.
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    % Score test using the five most common values of lambda.
    % In this case (given that all observations are positive the extended
    % score test for positive observations reduces to the standard score test
    % while that for negative is equal to NaN.
    [outSc]=ScoreYJpn(y,X);
%}


%{
    % Extended score test using Darwin data given by Yeo and Yohnson.
     y=[6.1, -8.4, 1.0, 2.0, 0.7, 2.9, 3.5, 5.1, 1.8, 3.6, 7.0, 3.0, 9.3, 7.5 -6.0]';
     n=length(y);
     X=ones(n,1);
     % Score and extended score test in the grid of lambda 1, 1.1, ..., 2
     la=[1:0.1:2];
     % Given that there are no explanatory variables the test must be
     % called with intercept 0
     outpn=ScoreYJpn(y,X,'intercept',0,'la',la);
     out=ScoreYJ(y,X,'intercept',0,'la',la);
     disp([la' outpn.Score(:,1) out.Score outpn.Score(:,2)])
%}



%% Beginning of code


nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

la=[-1 -0.5 0 0.5 1];

scoremle= false;
usefmin=true;

if nargin>2
    
    options=struct('la',la,'nocheck',0,'intercept',0,'scoremle',scoremle,...
        'usefmin',usefmin);
    
    UserOptions=varargin(1:2:length(varargin));
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:ScoreYJall:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    
    % Write in structure 'options' the options chosen by the user
    if nargin > 2
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
    end
    
    la=options.la;
    scoremle=options.scoremle;
    usefmin=options.usefmin;
    % SSElaMLE = residual sum of squares of the regression model based on
    % MLE of laPos and laNeg
    SSElaMLE=options.SSElaMLE;
end


% inizialize qunatities which do not depnd on la(i)
nonnegs = y >= 0;
negs = ~nonnegs;
ynonnegs=y(nonnegs);
ynegs=y(negs);

logG=sum(  sign(y) .* log(abs(y)+1)   )/n;
vneg=-ynegs+1;
vpos=ynonnegs+1;
logvpos=log(vpos);
logGpos=sum(logvpos)/n;
logvneg=log(vneg);
logGneg=sum(-logvneg)/ n;
G=exp(logG);
%Gpos=exp(logGpos);
%Gneg=exp(logGneg);
% Note that Gpos*Gneg=G


%  Sc= matrix lla-by-5 which contain:
% 1st col = global t test
% 2nd col = t test for pos
% 3rd col = t test for neg
% 4th col = F test for both
% 5th col = Lik. ratio test for both
%  for the values of \lambda specified in vector la
lla=length(la);
Sc=NaN(lla,5);
wini=NaN(n,1);


% The identity matrix of size p+1 can be
% computed once and for all
% eyep1=eye(p+1);

% loop over the values of \lambda
for i=1:lla
    z=y; % Initialized z and w
    wpos=wini;
    wneg=wini;
    w=wini;
    lai=la(i);
    Glaminus1=G^(lai-1);
    q=lai*Glaminus1;
    twomlambdai=2-lai;
    
    % Compute transformed values and constructed variables depending on lambda
    % transformation for non negative values
    if abs(lai)>1e-8  % if la is different from 0
        % vposlai=vpos.^lai;
        vposlai=exp(lai*logvpos);
        kpos= (1/lai+logGpos);
        znonnegs=(vposlai-1)/q;
        z(nonnegs)=znonnegs;
        wpos(nonnegs)=( vposlai  .*(logvpos-kpos)    +kpos) /q;
        wneg(nonnegs)=-znonnegs*logGneg;
        k= (1/lai+logG);
        w(nonnegs)=(  vposlai  .*(logvpos-k)    +k) /q;
    else % if la is equal to 0
        znonnegs=G*logvpos;
        z(nonnegs)=znonnegs;
        wpos(nonnegs)=znonnegs.*(logvpos/2 - logGpos );
        wneg(nonnegs)=-znonnegs*logGneg;
        w(nonnegs)=znonnegs.*(logvpos/2-logG);
    end
    
    % Transformation and constructed variables for negative values
    if   abs(lai-2)>1e-8 % la not equal 2
        % vnegtwomlambdai=vneg.^twomlambdai;
        vnegtwomlambdai=exp(twomlambdai*logvneg);
        qneg=twomlambdai* Glaminus1;
        kN=logGneg-1/twomlambdai;
        znegs=(1-vnegtwomlambdai )  /qneg;
        z(negs)=znegs;
        wpos(negs)=- znegs*logGpos;
        wneg(negs)=(vnegtwomlambdai .*(logvneg+kN) -kN)/qneg;
        
        k=logG-1/twomlambdai;
        w(negs)=(  vnegtwomlambdai .*  (logvneg+k)   -k  )/qneg;
        
    else  % la equals 2
        znegs=-logvneg/G;
        z(negs)=znegs;
        wpos(negs)=- znegs*logGpos;
        wneg(negs)=logvneg.*(logvneg/2 +logGneg)/ G;
        
        w(negs)=logvneg.*(logvneg/2+logG)/G;
    end
    
    % Compute residual sum of squares for null (reduced) model
    betaR=X\z;
    residualsR = z - X*betaR;
    % Sum of squares of residuals
    SSeR = norm(residualsR)^2;
    
    
    % Define augmented X matrix for overall constructed variable
    Xw=[X w];
    
    % New code
    beta=Xw\z;
    residuals = z - Xw*beta;
    % Sum of squares of residuals
    SSe = norm(residuals)^2;
    Ftestnum=(SSeR-SSe);
    Ftestden=SSe/(n-p-1);
    Ftest=Ftestnum/Ftestden;
    Sc(i,1)=-sign(beta(end))*sqrt(Ftest);
    
    
    
    % if vpos is empty all the observations are negative
    vposboo=~isempty(vpos);
    if vposboo == true
        % Define augmented X matrix
        Xw=[X wpos];
        
        % OLD slow code
        %         [Q,R] = qr(Xw,0);
        %         beta = R\(Q'*z);
        %         residuals = z - Xw*beta;
        %         % Sum of squares of residuals
        %         sse = norm(residuals)^2;
        %         % Compute t stat for constructed added variable
        %         ri = R\eyep1;
        %         xtxi = ri*ri';
        %         se = sqrt(diag(xtxi*sse/(n-p-1)));
        %         Sc(i,1) = -beta(end)/se(end);
        
        % New code
        beta=Xw\z;
        residuals = z - Xw*beta;
        % Sum of squares of residuals
        SSe = norm(residuals)^2;
        Ftestnum=(SSeR-SSe);
        Ftestden=SSe/(n-p-1);
        Ftest=Ftestnum/Ftestden;
        Sc(i,2)=-sign(beta(end))*sqrt(Ftest);
    else
        Sc(i,2) =NaN;
    end
    
    % if vneg is empty all the observations are negative
    vnegboo=~isempty(vneg);
    if vnegboo==true
        % Define augmented X matrix
        Xw=[X wneg];
        
        % OLD slow code
        %         [Q,R] = qr(Xw,0);
        %         beta = R\(Q'*z);
        %         residuals = z - Xw*beta;
        %         % Sum of squares of residuals
        %         sse = norm(residuals)^2;
        %         % Compute t stat for constructed added variable
        %         ri = R\eyep1;
        %         xtxi = ri*ri';
        %         se = sqrt(diag(xtxi*sse/(n-p-1)));
        %         Sc(i,2) = -beta(end)/se(end);
        
        % New code
        beta = Xw\z;
        residuals = z - Xw*beta;
        % Sum of squares of residuals
        SSe = norm(residuals)^2;
        Ftestnum=(SSeR-SSe);
        Ftestden=SSe/(n-p-1);
        Ftest=Ftestnum/Ftestden;
        Sc(i,3)=-sign(beta(end))*sqrt(Ftest);
    else
        Sc(i,3)=NaN;
    end
    
    % Compute the F test for the joint presence of both constructed
    % variables
    if vposboo==true && vnegboo==true
        % Compute residual sum of squares for full model
        Xww=[X wpos wneg];
        betaF=Xww\z;
        
        residualsF = z - Xww*betaF;
        % Sum of squares of residuals
        SSeF = norm(residualsF)^2;
        
        Ftestnum=(SSeR-SSeF)/2;
        Ftestden=SSeF/(n-p-2);
        Ftest=Ftestnum/Ftestden;
        Sc(i,4)=Ftest;
    elseif vposboo==true
        % If there are just positive observations F test is the square of
        % the t test for positive
        Sc(i,4)=Sc(i,2)^2;
    else % in this case there are just negative observations
        Sc(i,4)=Sc(i,3)^2;
    end
    
    
    if scoremle == true 
        if isempty(SSElaMLE)
            % compute the exact score test based on lik. ratio
            Likrat=ScoreYJmle(y,X,'la',lai,'sseReducedModel',SSeR,'usefmin',usefmin,'nocheck',1);
            Sc(i,5)=Likrat.Score;
            laMLE=Likrat.laMLE;
            SSElaMLE=Likrat.SSE;
        else
            Ftestnum=(SSeR-SSElaMLE)/2;
            Ftestden=SSElaMLE/(n-p-2);
            Ftest=Ftestnum/Ftestden;
            Sc(i,5)=Ftest;
        end
    end
end

% Store values of the score test inside structure outSC
outSC.Score=Sc;

% Also store MLE of lambda
if scoremle == true
    outSC.laMLE=laMLE;
end

end
%FScategory:REG-Transformations