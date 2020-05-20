function [outSC]=ScoreYJ(y,X,varargin)
%Computes the score test for Yeo and Johnson transformation
%
%<a href="matlab: docsearchFS('ScoreYJ')">Link to the help function</a>
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
%           la  :transformation parameter. Vector. It specifies for which values of the
%                 transformation parameter it is necessary to compute the
%                 score test.
%                 Default value of lambda is la=[-1 -0.5 0 0.5 1]; that
%                 is the five most common values of lambda
%               Example - 'la',[0 0.5]
%               Data Types - double
%
%           Lik : likelihood for the augmented model. Scalar.
%                   If 1 the value of the likelihood for the augmented model will be produced
%                 else (default) only the value of the score test will be given
%               Example - 'Lik',0
%               Data Types - double
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
%        outSC.Score    =    score test. Scalar. t test for additional
%                            constructed variable
%        outSC.Lik      =    value of the likelihood. Scalar. This output
%                           is produced only if input value Lik =1
%
% See also: FSRfan, Score, normBoxCox, normYJ, ScoreYJpn
%
% References:
%
% Yeo, I.K. and Johnson, R. (2000), A new family of power
% transformations to improve normality or symmetry, "Biometrika", Vol. 87,
% pp. 954-959.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('ScoreYJ')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples


%{
    % Score with all default options for the wool data.
    % Load the wool data.
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    % Score test using the five most common values of lambda
    [outSc]=ScoreYJ(y,X);
    disp(outSc.Score)
%}

%{
    % Score with optional arguments.
    % Loyalty cards data.
    load('loyalty.txt');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    % la = vector containing the values of the transformation
    % parameters which have to be tested
    la=[0.25 1/3 0.4 0.5];
    [outSc]=ScoreYJ(y,X,'la',la,'intercept',1);
%}

%{
    % Compare Score test using BoxCox and YeoJohnson for the wool data.
    % Wool data.
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    % Define vector of transformation parameters
    la=[-1:0.01:1];
    % Score test using YeoJohnson transformation
    [outYJ]=ScoreYJ(y,X,'la',la);
    % Score test using YeoJohnson transformation
    [outBC]=Score(y,X,'la',la);
    plot(la',[outBC.Score outYJ.Score])
    xlabel('\lambda')
    ylabel('Value of the score test')
    legend({'BoxCox' 'YeoJohnson'})
%}

%{
    %% Score test using Darwin data given by Yeo and Yohnson.
     y=[6.1, -8.4, 1.0, 2.0, 0.7, 2.9, 3.5, 5.1, 1.8, 3.6, 7.0, 3.0, 9.3, 7.5 -6.0]';
     n=length(y);
     X=ones(n,1);
     la=-1:0.01:2.5;
     % Given that there are no explanatory variables the test must be
     % called with intercept 0
     out=ScoreYJ(y,X,'intercept',0,'la',la,'Lik',1);
     plot(la',out.Score)
     xax=axis;
     line(xax(1:2),zeros(1,2))
     xlabel('\lambda')
     ylabel('Value of the score test')
     title('Value of the score test is 0 in correspondence of $\hat \lambda =1.305$','Interpreter','Latex')
%}

%% Beginning of code
nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

la=[-1 -0.5 0 0.5 1];
Likopt=0;

if nargin>2
    
    options=struct('Lik',0,'la',la,'nocheck',0,'intercept',0); % ,'mingreater0',mingreater0);
    
    UserOptions=varargin(1:2:length(varargin));
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:Score:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
    % mingreater0=options.mingreater0;
    Likopt=options.Lik;
end


%  Sc= vector which contains the t test for constructed variables for the
%  values of \lambda specified in vector la
lla=length(la);
Sc=zeros(lla,1);

% Lik is a vector which contains the likelihoods for diff. values of la
if Likopt==1
    Lik=Sc;
end

% value related to the Jacobian
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
    z=y; % Initialized z and w
    w=y;
    
    lai=la(i);
    % Glaminus1=G^(lai-1);
    Glaminus1=exp((lai-1)*logG);
    
    % Define transformed and constructed variable
    % transformation for non negative values
    % Compute transformed values and constructed variables depending on lambda
    % transformation for non negative values
    if abs(lai)>1e-8  % if la is different from 0
        q=lai*Glaminus1;
        % vposlai=vpos.^lai;
        vposlai=exp(lai*logvpos);
        znonnegs=(vposlai-1)/q;
        z(nonnegs)=znonnegs;
        k= (1/lai+logG);
        w(nonnegs)=(  vposlai  .*(logvpos-k)    +k) /q;
    else % if la is equal to 0
        znonnegs=G*logvpos;
        z(nonnegs)=znonnegs;
        w(nonnegs)=znonnegs.*(logvpos/2-logG);
    end
    
    % Transformation and constructed variables for negative values
    if   abs(lai-2)>1e-8 % la not equal 2
        twomlambdai=2-lai;
        % vnegtwomlambdai=vneg.^twomlambdai;
        vnegtwomlambdai=exp(twomlambdai*logvneg);
        qneg=twomlambdai* Glaminus1;
        znegs=(1-vnegtwomlambdai )  /qneg;
        z(negs)=znegs;
        
        k=logG-1/twomlambdai;
        w(negs)=(  vnegtwomlambdai .*  (logvneg+k)   -k  )/qneg;
        
    else  % la equals 2
        znegs=-logvneg/G;
        z(negs)=znegs;
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
    Sc(i)=-sign(beta(end))*sqrt(Ftest);
    
    % Store the value of the likelihood for the model which also contains
    % the constructed variable
    if Likopt==1
        Lik(i)=-n*log(SSe/n);
    end
end

% Store values of the score test inside structure outSC
outSC.Score=Sc;

% Store values of the likelihood inside structure outSC
if Likopt==1
    outSC.Lik=Lik;
end

end
%FScategory:REG-Transformations