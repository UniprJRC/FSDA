function [out] = simulateLM(n,varargin)
%simulateLM simulates linear regression data with prespecified values of statistical indexes.
%
%<a href="matlab: docsearchFS('simulateLM')">Link to the help function</a>
%
% simulateLM simulates linear regression data. It is possible to specify:
% 1) the requested value of R2 (or equivaletly its SNR);
% 2) the values of the beta coefficients (possibly sparse);
% 3) the correlation (covariance) matrix among the explanatory variables.
% 4) the value of the intercept term.
% 5) the distribution to use to generate the Xs;
% 6) the distribution to use to generate the ys.
% 7) the MSOM contamination in Xs and ys.
% 8) the VIOM contamination in ys.
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
%
%     SNR  : Signal to noise ratio characterizing the simulation. This
%            is defined such that sigma_error == sqrt(var(X_u*beta_true)/SNR)
%            The default is SNR=='' and R2 is used instead.
%            Example - 'SNR',10
%            Data Types - double
%
%     beta :   the values of the beta coefficients. Vector. Vector which
%              contains the values of the regression coefficients. The
%              default is a vector of ones.
%                 Example - 'beta',[3 5 8]
%                 Data Types - double
%
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
%       nexpl   : number of explanatory variables. If vector beta is
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
%       pMSOM   : Proportion of MSOM outliers. The default is 10% MSOM 
%                 contmaination.
%                 Example - 'pMSOM',0.25
%                 Data Types - double
% 
%       pVIOM   : Proportion of VIOM outliers (non-overlapping with MSOM). 
%                 The default is 10% VIOM contmaination.
%                 Example - 'pVIOM',0.25
%                 Data Types - double
% 
%   shiftMSOMe  : Mean-shift on the error terms for MSOM outliers.
%                 Default value shiftMSOMe==10.
%                 Example - 'shiftMSOMe',-3
%                 Data Types - double
% 
%   predxMSOM   : Predictors subject to a mean shift by MSOM. It is a 
%                 p-dimensional vector indexing design matrix columns.
%                 Default value is to contaminate only the non-zero
%                 entries of beta_true (excluding the intercept).
%                 Example - 'predMSOM',true(2,1)
%                 Data Types - boolean
% 
%   shiftMSOMx  : Mean-shift on the predictor terms for MSOM outliers.
%                 Default value shiftMSOMx==10.
%                 Example - 'shiftMSOMx',3
%                 Data Types - double
% 
%   inflVIOMe   : Variance-inflation for the errors subject to a VIOM.
%                 Default value is inflVIOMe==10.
%                 Example - 'inflVIOMe',5
%                 Data Types - double
%
%
% Output:
%
%         out:   structure which contains the following fields
%
%    out.y        = simulated response. Vector. Column vector of length n
%                   containing the response.
%    out.X        = simulated regressors. Matrix . Matrix of size
%                   n-times-nexpl containing the values of the regressors.
%
%  Optional Output (for pVIOM+pMSOM>0):
%
%   out.yc        = Contaminated response vector.
%   out.Xc        = Contaminated response vector.
% out.ind_clean   = Indexes for non-outlying cases.
% out.ind_MSOM    = Indexes for MSOM outlying cases.
% out.ind_VIOM    = Indexes for VIOM outlying cases.
% out.vareps      = Variance for the uncontaminated errors.
%
%
% See also simulateTS
%
% References:
%
% Insolia, L., F. Chiaromonte, and M. Riani (2020a). 
% “A Robust Estimation Approach for Mean-Shift and Variance-Inflation Outliers”. 
% In press.
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
    out=simulateLM(100,'plots',true);
%}

%{
    %% Simulate with prefixed value of R2.
    % Set value of R2;
    R2=0.82;
    n=10000;
    out=simulateLM(n,'R2',R2);
    outLM=fitlm(out.X,out.y);
%}

%{
    %% Use prefixed correlation matrix for cov(X).
    % Set value of R2;
    R2=0.26;
    n=10000;
    A = gallery('moler',5,0.2);
    out=simulateLM(n,'R2',R2,'SigmaX',A);
    outLM=fitlm(out.X,out.y)
%}

%{

    %% Use prefixed values of R2, beta and intercept.
    % Set value of R2.
    R2=0.92;
    beta=[3; 4; 5; 2; 7];
    intercept=43;
    n=100000;
    out=simulateLM(n,'R2',R2,'beta',beta);
    outLM=fitlm(out.X,out.y);

%}

%{
    % Sim study.
    % Compare the distribution of values of R2 with data generated from 
    % Normal with those generated from Student T with 5 degrees of freedom.
    % Set value of R2.
    R2=0.92;
    beta=[3; 4; 5; 2; 7; 2; 3];
    nsimul=1000;
    R2all=zeros(nsimul,2);
    n=100;
    df=5;
    for j=1:nsimul
        % Data generated from Normal
        out=simulateLM(n,'R2',R2,'beta',beta);
        outLM=fitlm(out.X,out.y);
        R2all(j,1)=outLM.Rsquared.Ordinary;
        % Data generated from T(5)
        out=simulateLM(n,'R2',R2,'beta',beta,'distriby','T','distribypars',df);
        outLM=fitlm(out.X,out.y);
        R2all(j,2)=outLM.Rsquared.Ordinary;
    end
    boxplot(R2all,'Labels',{'Normal', 'T(5)'});
%}

%{

    %% Use SNR and include MSOM (on active features) and VIOM contamination
    SNR=3;
    beta=[2, 2, 0, 0];
    intercept=1;
    n=100;
    out=simulateLM(n,'SNR',SNR,'beta',beta, 'pMSOM', 0.1, 'pVIOM', 0.2, 'plots', 1);
    X = out.X;
    y = out.y;
    outLM=fitlm(X,y);
    Xc = out.Xc;
    yc = out.yc;
    outLM2=fitlm(Xc,yc);
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
SNR = '';
pMSOM = 0;
pVIOM = 0;
shiftMSOMe = -10;
predxMSOM = '';
shiftMSOMx = 10;
inflVIOMe = 10;


options=struct('R2',R2,...
    'beta',beta,'SigmaX',SigmaX,...
    'distribX',distribX,'distribXpars',distribXpars,...
    'distriby',distriby,'distribypars',distribypars,...
    'nexpl',nexpl,'intercept',intercept,'plots',plots, ...
    'SNR', SNR, 'pMSOM', pMSOM, 'pVIOM', pVIOM, ...
    'shiftMSOMe', shiftMSOMe, 'predxMSOM', predxMSOM, ...
    'shiftMSOMx', shiftMSOMx, 'inflVIOMe', inflVIOMe);


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
    SNR = options.SNR;
    pMSOM = options.pMSOM;
    pVIOM = options.pVIOM;
    shiftMSOMe = options.shiftMSOMe;
    if isempty(predxMSOM)
        predxMSOM= abs(beta)>0;
    else
        predxMSOM = options.predxMSOM;
    end
    shiftMSOMx = options.shiftMSOMx;
    inflVIOMe = options.inflVIOMe;

    if R2<0 || R2>=1
        error('FSDA:simulateLM:WrongOpt','R2 must belong to [0, 1)');
    end
    if SNR<=0
        error('FSDA:simulateLM:WrongOpt','SNR must be greater than zero');
    end
    
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
    
    
    if ischar(distribX)
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
    else
        % In this case the user has directly supplied matrix X.
        % Make sure that the size of X is n-by-nexpl
        X=distribX;
        [nchk,nexplchk]=size(X);
         if nchk~=n
            error('FSDA:simulateLM:WrongOpt',['supplied matrix X must have  '  ...
                num2str(n) ' rows']);
         end
         if nexpl~=nexplchk
            error('FSDA:simulateLM:WrongOpt',['supplied matrix X must have '  ...
                num2str(nexpl) ' columns']);
         end
    end
    
    
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
    
    p=nexpl+intercept;
    
    % Divide by std and multyply by a small sample correction factor.
     err=sqrt((n)/(n-p))*err/std(err,1);
    % err=err/std(err,1);
    
    if R2>0 && isempty(SNR)
        % Find var(\epsilon) which produces a value of R2 centered around
        % the one which has been requested.
        vareps=(intercept+beta'*SigmaX*beta)*((1 - R2)/R2);
        y=intercept+X*beta(:)+err*sqrt(vareps);
    elseif ~isempty(SNR)
        vareps = var(X*beta(:)) / SNR;
        y=intercept+X*beta(:)+err*sqrt(vareps);
    else
        y=intercept+err;
    end
    
    if plots==true
        yXplot(y,X);
    end
end

%% checks on the explanatory variables


%% add contamination as MSOM and VIOM
    
if pMSOM + pVIOM >0
    
    eu = err*sqrt(vareps);
    Xu = X;

    % number of units arising from a MSOM
    nMSOM = floor(n * pMSOM);
    % indexes MSOM
    indMSOM = randperm(n, nMSOM);
    % MSOM (also on X)
    Xc = Xu;
    Xc(indMSOM, predxMSOM) = Xu(indMSOM, predxMSOM) + shiftMSOMx; 
    ec = eu;
    ec(indMSOM) = eu(indMSOM) + shiftMSOMe;

    % number of units arising from a VIOM
    nVIOM = floor(n * pVIOM);
    % indexes VIOM
    indVIOM = datasample(setdiff(1:n, indMSOM), nVIOM, 'replace', false);
    %VIOM
    ec(indVIOM) = eu(indVIOM) .* sqrt(inflVIOMe);
    % contaminated model
    yc = Xu * beta(:) + ec;

    % total contaminated units
    indcont = unique([indMSOM, indVIOM]);
    % clean obs indexes
    indkeep = setdiff(1:n, indcont);
    
    if plots==true
        indunit = zeros(n, 1);
        indunit(indMSOM) = 1;
        indunit(indVIOM) = 2;
        yXplot(yc,Xc, indunit);
    end

    %% save output

    out.Xc = Xc;
    out.yc = yc;
    out.ind_clean = indkeep;
    out.ind_MSOM = indMSOM;
    out.ind_VIOM = indVIOM;
    out.vareps = vareps;
    
end

out.X = X;
out.y = y;

end

%FScategory:REG-Regression