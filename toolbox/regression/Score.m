function [outSC]=Score(y,X,varargin)
%Score computes the score test for Box-Cox transformation
%
%<a href="matlab: docsearchFS('Score')">Link to the help function</a>
%
%  Required input arguments:
%
%    y:         Response variable. Vector. A vector with n elements that
%               contains the response variable.
%               It can be either a row or a column vector.
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
%           la : transformation parameter. Vector. It specifies for which values of the
%                 transformation parameter it is necessary to compute the
%                 score test.
%                 Default value of lambda is la=[-1 -0.5 0 0.5 1]; that
%                 is the five most common values of lambda.
%               Example - 'la',[0 0.5]
%               Data Types - double
%
%           Lik : likelihood for the augmented model. Boolean.
%                   If true the value of the likelihood for the augmented
%                   model will be produced,
%                 else (default) only the value of the score test will be
%                 given.
%               Example - 'Lik',false
%               Data Types - logical
%
%       nocheck : Check input arguments. Boolean.
%               If nocheck is equal to true no check is performed on
%                 matrix y and matrix X. Notice that y and X are left
%                 unchanged. In other words, the additional column of ones
%                 for the intercept is not added. As default nocheck=false.
%               Example - 'nocheck',true
%               Data Types - boolean
%
%      tukey1df : Tukey's one df test. Boolean.
%                 Tukey's one degree of freedome test for non-additivity.
%                 The constructed variable is given by
%                 \[
%                  w_T(\lambda)= (\hat z(\lambda) - \overline  z(\lambda))^2 / 2 \overline  z(\lambda)
%                 \]
%                 where $z(\lambda)$ is the transformed response, and
%                 $\hat z(\lambda)$ are the fitted values on the
%                 transformed response. The t test on the constructed
%                 variable above provides a test from departures from the
%                 assumed linear model and is known in the literature as
%                 Tukey's one degree of freedome test for non-additivity.
%                 If tukey1df is true the test is computed and returned
%                 inside output structure with fieldname ScoreT else
%                 (default) the value of the test is not computed.
%               Example - 'tukey1df',true
%               Data Types - boolean
%
%  Output:
%
%  The output consists of a structure 'outSC' containing the following fields:
%
%        outSC.Score    =    score test. Vector. Vector of length
%                            length(lambda) which contains the value of the
%                            score test for each value of lambda specified
%                            in optional input parameter la. If la is not
%                            specified, the vector will be of length 5 and
%                            contains the values of the score test for the
%                            5 most common values of lambda.
%        outSC.ScoreT   =    value of the Tukey's one degree of freedome
%                            test for non-additivity. Scalar. This output
%                            is produced only if optional input tukey1df is
%                            true.
%        outSC.Lik      =    value of the likelihood. Scalar. This output
%                           is produced only if optional input Lik=1.
%
% See also: FSRfan, ScoreYJ, ScoreYJpn, normBoxCox, normYJ
%
% References:
%
% Atkinson, A.C. and Riani, M. (2000), "Robust Diagnostic Regression
% Analysis", Springer Verlag, New York. [see equation 2.30 for the
% expression for score test statistic]
%
% Copyright 2008-2024.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('Score')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples

%{
    %% Score with all default options.
    % Wool data.
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    % Score test using the five most common values of lambda.
    [outSC]=Score(y,X);
    disp('Values of the score test')
    disp({'la=-1' 'la=-0.5' 'la=0' 'la=0.5' 'la=1'})
    disp(outSC.Score')
%}

%{
    % Score with optional arguments.
    % Loyalty cards data.
    load('loyalty.txt');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    % la = vector containing the values of the transformation
    % parameters that have to be tested.
    la=[0.25 1/3 0.4 0.5];
    [outSc]=Score(y,X,'la',la,'intercept',false);
%}




%% Beginning of code

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = aux.chkinputR(y,X,nnargin,vvarargin);

% Add the extra check on vector y
if min(y)<0
    error('FSDA:Score:ynegative','Score test using BoxCox cannot be computed because min(y) is smaller than 0. Please use Yeo-Johnson family')
end

Likboo=false;
tukey1df=false;
la=[-1 -0.5 0 0.5 1];
if coder.target('MATLAB')

    options=struct('Lik',Likboo,'la',la,'nocheck',false,'intercept',false,'tukey1df',tukey1df);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid.
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:Score:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options.
        aux.chkoptions(options,UserOptions)
    end
end

% Write in structure 'options' the options chosen by the user.
if nargin > 2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    la=options.la;
    Likboo=options.Lik;
    tukey1df=options.tukey1df;
end


%  Sc= vector which contains the t test for constructed variables for the
%  values of \lambda specified in vector la.
lla=length(la);
Sc=zeros(lla,1);
if tukey1df == true
    ScT=Sc;
end

% Lik is a vector which contains the likelihoods for diff. values of la.
Lik=Sc;

% Geometric mean of the observations.
logy=log(y);
G=exp(sum(logy)/n);
logG=log(G);

% loop over the values of \lambda.
for i=1:lla
    lai=la(i);
    % Define transformed and constructed variable
    if abs(lai)<1e-8
        z=G*logy;
        w=G*logy.*(logy/2-logG);
    else
        % laiGlaim1=lai*G^(lai-1);
        laiGlaim1 =lai*exp((lai-1)*logG);
        % ylai=y.^lai;
        ylai=exp(lai*logy);
        ylaim1=ylai-1;
        z=ylaim1/laiGlaim1;
        w=(ylai.*logy-ylaim1*(1/lai+logG))/laiGlaim1;

        % OLD slow code
        % z=(y.^la(i)-1)/(la(i)*G^(la(i)-1));
        % w=(y.^la(i).*log(y)-(y.^la(i)-1)*(1/la(i)+log(G)))/(la(i)*G^(la(i)-1));
    end

    % Compute tstat on costructed variable
    eyepplus1=eye(p+1);
    [tstatw,sse]=ScoreCore(z,X,w,n,p,eyepplus1);
    Sc(i)=tstatw;
    if tukey1df == true
        b=X\z;
        zhat=X*b;
        mz=sum(z)/n;
        w=(zhat -mz).^2/(2*mz);
        tstatw=ScoreCore(z,X,w,n,p,eyepplus1);
        ScT(i)=tstatw; 
    end

    % Store the value of the likelihood for the model which also contains
    % the constructed variable.
    if Likboo==true
        Lik(i)=-n*log(sse/n);
    end
end

% Store values of the score test inside structure outSC.
outSC.Score=Sc;
if tukey1df == true
    outSC.ScoreT=ScT;
end

% Store values of the likelihood inside structure outSC
if Likboo==true
    outSC.Lik=Lik;
else
    outSC.Lik=NaN;
end

end

function [tstatw,sse]=ScoreCore(z,X,w,n,p,eyepplus1)
% Define augmented X matrix.
Xw=[X w];
[Q,R] = qr(Xw,0);
beta = R\(Q'*z);
residuals = z - Xw*beta;
% Sum of squares of residuals.
sse = norm(residuals)^2;
% Compute t stat of constructed added variable.
ri = R\eyepplus1;
xtxi = ri*ri';
se = sqrt(diag(xtxi*sse/(n-p-1)));
tstatw= -beta(end)/se(end);
end

%FScategory:REG-Transformations