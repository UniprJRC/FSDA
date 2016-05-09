function [outSC]=Score(y,X,varargin)
%Score computes the score test for transformation
%
%<a href="matlab: docsearchFS('Score')">Link to the help function</a>
%
%  Required input arguments:
%
%    y: A vector with n elements that contains the response variable. y can
%       be both a row of column vector.
%    X: Data matrix of explanatory variables (also called 'regressors') of
%       dimension (n x p-1). Rows of X represent observations, and columns
%       represent variables.
%       Missing values (NaN's) and infinite values (Inf's) are allowed,
%       since observations (rows) with missing or infinite values will
%       automatically be excluded from the computations.
%
%  Optional input arguments:
%
%  intercept :  Indicator for constant term. Scalar. If 1, a model with
%               constant term will be fitted (default), if 0, no constant
%               term will be included.
%               Example - 'intercept',1 
%               Data Types - double
%           la  :It specifies for which values of the
%                 transformation parameter it is necessary to compute the
%                 score test. Vector.
%                 Default value of lambda is la=[-1 -0.5 0 0.5 1]; that
%                 is the five most common values of lambda
%               Example - 'la',[0 0.5]
%               Data Types - double
%           Lik : likelihood for the augmented model. Scalar.
%                   If 1 the value of the likelihood for the augmented model will be produced
%                 else (default) only the value of the score test will be given
%               Example - 'Lik',0
%               Data Types - double
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
%                           is produced just if input value Lik =1
%
% See also
% 
% FSRfan
%
% References:
%
% Atkinson Riani (2000), equation (2.30) for the expression
% for score test statistic.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('Score')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples

%{
    % Score with all default options.
    % Wool data.
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    % Score test using the five most common values of lambda
    [outSc]=Score(y,X);
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
    [outSc]=Score(y,X,'la',la,'intercept',0);
%}

%{

%}


nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% Add the extra check on vector y
if min(y)<0
    error('FSDA:Score:ynegative','Score test cannot be computed because min(y) is smaller than 0')
end

options=struct('Lik',0,'la',[-1 -0.5 0 0.5 1],'nocheck',0,'intercept',0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:Score:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end

la=options.la;


%  Sc= vector which contains the t test for constructed variables for the
%  values of \lambda specified in vector la
lla=length(la);
Sc=zeros(lla,1);

% Lik is a vector which contains the likelihoods for diff. values of la
Lik=Sc;

% Geometric mean of the observations
G=exp(mean(log(y)));

% loop over the values of \lambda
for i=1:lla;
  
    % Define transformed and constructed variable
    if abs(la(i))<1e-8;
        z=G*log(y);
        w=G*log(y).*(log(y)/2-log(G));
    else
        z=(y.^la(i)-1)/(la(i)*G^(la(i)-1));
        w=(y.^la(i).*log(y)-(y.^la(i)-1)*(1/la(i)+log(G)))/(la(i)*G^(la(i)-1));
    end
    
    % Define augmented X matrix
    Xw=[X w];
    [Q,R] = qr(Xw,0);
    beta = R\(Q'*z);
    residuals = z - Xw*beta;
    % Sum of squares of residuals
    sse = norm(residuals)^2;
    % Compute t stat for constructed added variable
    ri = R\eye(p+1);
    xtxi = ri*ri';
    se = sqrt(diag(xtxi*sse/(n-p-1)));
    Sc(i) = -beta(end)/se(end);
    
    % Store the value of the likelihood for the model which also contains
    % the constructed variable
    if options.Lik==1
        Lik(i)=-n*log(sse/n);
    end
end

% Store values of the score test inside structure outSC
outSC.Score=Sc;

% Store values of the likelihood inside structure outSC
if options.Lik==1
    outSC.Lik=Lik;
end

end
%FScategory:REG-Transformations