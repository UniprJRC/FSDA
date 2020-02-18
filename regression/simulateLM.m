function [y, X] = simulateLM(n,varargin)
%Simulate linear regression data
%
%<a href="matlab: docsearchFS('simulateLM')">Link to the help function</a>
%
% simulateLM simulates linear regression data. It is possible to specify:
% 1) the requested value of R2;
% 2) the values of the beta coefficients;
% 3) the correlation (covariance) matrix among the explanatory variables.
% 4) the value of the intercept term.
% 5) the distribution to use to generate the Xs;
% 6) the distribution to use to generate the ys.
%
%  Required input arguments:
%
%         n  :  sample size. Scalar. n is a positive integer
%               which defines the length of the simulated data. For example
%               if n=100, y will be 100x1 and X will be 100xp.
%
%  Optional input arguments:
%
%      R2 :  Squared multiple correlation coefficient (R2). Scalar. The
%            requested value of R2. A number in the
%            interval [0 1] which specifies the requested value of R2.
%            The default is to simulate regression data with R2=0;
%                 Example - 'R2',0.90
%                 Data Types - double
%     beta :   the values of the beta coefficients. Vector. Vector which
%              contains the values of the regression coefficients. The
%              default is a vector of ones.
%                 Example - 'beta',[3 5 8]
%                 Data Types - double
%    SigmaX :   the correlation matrix. Matrix. Positive definite matrix
%               which contains the correlation matrix among regressors. The
%               default is the identity matrix.
%                 Example - 'Sigma', gallery('lehmer',5)
%                 Data Types - double
%   distribX : distribution to use to simulate the regressors. Character.
%              Character which specifies the distribution to use to
%              simulate the values of the explanatory variables.
%              For the list of valid names see MATLAB function random.
%              Default is to use the Standard normal distribution.
%                 Example - 'distribX', 'Beta'
%                 Data Types - double
% distribXpars : parameters of the distribution to use in distribX. Vector.
%              Scalar value or array of scalar values containing the
%              distribution parameters specified in distribX.
%                 Example - 'distribXpars', '[0.2 0.6]'
%                 Data Types - double
%   distriby : distribution to use to simulate the response. Character.
%              Character which specifies the distribution to use to
%              simulate the values of the explanatory variables. The
%              default is to use the Standard normal distribution.
%                 Example - 'distriby', 'Lognormal'
%                 Data Types - double
% distribypars : parameters of the distribution to use in distriby. Vector.
%              Scalar value or array of scalar values containing the
%              distribution parameters specified in distriby. For examples
%              if distriby is 'Lognormal' and 'distribypars' is [2 10], the
%              errors are generated according to a Log Normal distribution
%              with parameters mu and sigma respectively equal to 2 and 10.
%                 Example - 'distribypars', '[2 10]'
%                 Data Types - double
%       nexpl   : number of explanatory varibles. If vector beta is
%                 supplied nexpl is equal to length(beta). Similarly if
%                 sigmaX is supplied nexpl is set equal to size(sigmaX,1).
%                 Note that both nexpl is supplied together with beta and SigmaX it is check that
%                 nexpl =length(beta) = size(SigmaX,1). If options beta and
%                 sigmaX are empty nexpl is set equal to 3.
%                 Example - 'distribypars', '[2 10]'
%                 Data Types - double
%    intercept : value of the intercept to use. Scalar. The default value
%               for intercept is 0.
%                 Example - 'intercept', '10'
%                 Data Types - double
%       plots : Plot on the screen. Boolean.
%               If plots = true, the yXplot which shows the response
%               against all the explanatory variables s shown on the
%               screen. The default value for plots is false, that is no
%               plot is shown on the screen.
%                 Example - 'plots',false
%                 Data Types - single | double
%
%
%  Output:
%
%           y :  simulated response. Vector. Column vector of length n
%               containing the response.
%           X :  simulated regressors. Matrix . Matrix of size
%                n-times-nexpl containing the values of the regressors.
%
% See also simulateTS
%
% References:
%
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('simulateLM')">Link to the help function</a>
%
%$LastChangedDate:: 2019-05-14 16:04:25 #$: Date of the last commit


% Examples:

%{
    %% Use all defaul options.
    % Simulate 100 observations y and X (uncorrelated with y) using standard normal distribution.
    [y,X]=simulateLM(100,'plots',true);
%}

%{
    %% Simulate with prefixed value of R2.
    % Set value of R2;
    R2=0.82;
    n=10000;
    [y,X]=simulateLM(n,'R2',R2);
    out=fitlm(X,y);
%}

%{
    %% Use prefixed correlation matrix for cov(X).
    % Set value of R2;
    R2=0.26;
    n=10000;
    A = gallery('moler',5,0.2);
    [y,X]=simulateLM(n,'R2',R2,'SigmaX',A);
    out=fitlm(X,y)
%}

%{
    %% Use prefixed values of R2, beta and intercept.
    % Set value of R2.
    R2=0.92;
    beta=[3; 4; 5; 2; 7];
    intercept=43;
    n=100000;
    [y,X]=simulateLM(n,'R2',R2,'beta',beta);
    out=fitlm(X,y);
%}

%{
    % Sim study.
    % Compare the distribution of values of R2 with data generated from 
    % Normal with those generated from Student T with 5 degrees of freedom.
    nsimul=1000;
    R2all=zeros(nsimul,2);
    n=100;
    df=5;
    for j=1:nsimul
        % Data generated from Normal
        [y,X]=simulateLM(n,'R2',R2,'beta',beta);
        out=fitlm(X,y);
        R2all(j,1)=out.Rsquared.Ordinary;
        % Data generated from T(5)
        [y,X]=simulateLM(n,'R2',R2,'beta',beta,'distriby','T','distribypars',df);
        out=fitlm(X,y);
        R2all(j,2)=out.Rsquared.Ordinary;
    end
    boxplot(R2all,'Labels',{'Normal', 'T(5)'});
%}

%% Beginning of code
if nargin<1
    error('FSDA:simulateLM:MissingInputs','Input number of observations is missing');
end

R2=0;
p=3;
nexpl=p;
beta=ones(p,1);
SigmaX=eye(p);
distribX = 'normal';
distribXpars=[0 1];
distriby = 'normal';
distribypars = [0 1];
plots=false;
intercept=0;


options=struct('R2',R2,...
    'beta',beta,'SigmaX',SigmaX,...
    'distribX',distribX,'distribXpars',distribXpars,...
    'distriby',distriby,'distribypars',distribypars,...
    'nexpl',nexpl,'intercept',intercept,'plots',plots);


%% User options


UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:simulateLM:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present in
    % structure options Remark: the nocheck option has already been dealt
    % by routine chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:simulateLM:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
    
    % Check the presence of input options beta, SigmaX and nexpl
    betaboo=max(strcmp(UserOptions,'beta'))==1;
    SigmaXboo=max(strcmp(UserOptions,'SigmaX'))==1;
    nexplboo=max(strcmp(UserOptions,'nexpl'))==1;
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    R2=options.R2;
    nexpl=options.nexpl;
    beta=options.beta;
    SigmaX=options.SigmaX;
    distribX = options.distribX;
    distribXpars=options.distribXpars;
    distriby = options.distriby;
    distribypars = options.distribypars;
    plots=options.plots;
    intercept=options.intercept;
    
    
    % Preliminary checks both beta, nexpl and sigmaX have been supplied
    if betaboo==true && SigmaXboo==true && nexplboo==true
        
        if nexpl~=size(SigmaX,1)
            error('FSDA:simulateLM:WrongOpt',['Length of supplied vector beta ' ...
                'must be equal to number of rows (columns) of matrix SigmaX']);
        end
        
        if nexpl~=length(beta)
            error('FSDA:simulateLM:WrongOpt',['Length of supplied vector beta ' ...
                'must be equal to number of rows (columns) of matrix SigmaX']);
        end
    end
    
    % Preliminary checks just beta and sigmaX have been supplied
    if betaboo==true && SigmaXboo==true && nexplboo==false
        nexpl=length(betaboo);
        if nexpl~=size(SigmaXboo,1)
            error('FSDA:simulateLM:WrongOpt',['Length of supplied vector beta ' ...
                'must be equal to number of rows (columns) of matrix SigmaX']);
        end
    end
    
    % Preliminary checks just beta and nexpl have been supplied
    if betaboo==true && SigmaXboo==false  && nexpl==true
        nexpl=length(beta);
        if nexpl~=length(beta)
            error('FSDA:simulateLM:WrongOpt',['Length of supplied vector beta ' ...
                'must be equal to input option nexpl']);
        end
        SigmaX=eye(nexpl);
    end
    
    % Preliminary checks just SigmaX and nexpl have been supplied
    if betaboo==false && SigmaXboo==true && nexplboo==true
        nexplchk=size(SigmaXboo,1);
        if nexpl~=nexplchk
            error('FSDA:simulateLM:WrongOpt',['nexpl ' ...
                'must be equal to number of rows (columns) of matrix SigmaX']);
        end
        beta=ones(nexpl,1);
    end
    
    % Preliminary checks just beta has been supplied
    if betaboo==true &&  SigmaXboo == false && nexplboo==false
        nexpl=length(beta);
        SigmaX=eye(nexpl);
    end
    
    % Preliminary checks just SigmaX has been supplied
    if betaboo==false &&  SigmaXboo == true && nexplboo==false
        nexpl=size(SigmaX,1);
        beta=ones(nexpl,1);
    end
    
    % Preliminary checks just nexpl has been supplied
    if betaboo==false &&  SigmaXboo == false && nexplboo==true
        beta=ones(nexpl,1);
        SigmaX=eye(nexpl);
    end
    
    [T,err] = cholcov(SigmaX);
    if err ~= 0
        error('FSDA:mvnrnd:BadCovariance2DSymPos','WrongSigma');
    end
    lXpars=length(distribXpars);
    if lXpars==1
        X = random(distribX,distribXpars,n,nexpl);
    elseif lXpars==2
        X = random(distribX,distribXpars(1),distribXpars(2),n,nexpl);
    elseif lXpars==3
        X = random(distribX,distribXpars(1),distribXpars(2),distribXpars(3),n,nexpl);
    else
        X = random(distribX,distribXpars(1),distribXpars(2),distribXpars(3),distribXpars(4),n,nexpl);
    end
    % Generate the X in such a way their corr is SigmaX
    X=X*T;
    
    lypars=length(distribypars);
    if lypars==1
        err = random(distriby,distribypars,n,1);
    elseif lypars==2
        err = random(distriby,distribypars(1),distribypars(2),n,1);
    elseif lypars==3
        err = random(distriby,distribypars(1),distribypars(2),distribypars(3),n,1);
    else
        err = random(distriby,distribypars(1),distribypars(2),distribypars(3),distribypars(4),n,1);
    end
    
    % Divide by std and multyply by a small sample correction factor.
    err=sqrt((n)/(n-p))*err/std(err,1);
    
    if R2>0
        % Find var(\epsilon) which produces a value of R2 centered around
        % the one which has been requested.
        vareps=(beta'*SigmaX*beta)*((1 - R2)/R2);
        
        y=intercept+X*beta(:)+err*sqrt(vareps);
    else
        y=intercept+err;
    end
    
    if plots==true
        yXplot(y,X);
    end
end

%% checks on the explanatory variables

end

%FScategory:REG-Regression