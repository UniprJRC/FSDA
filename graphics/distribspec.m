function [p, h] = distribspec(pd, specs, region, varargin)
%distribspec plots a probability density function between specification limits
%
% <a href="matlab: docsearchFS('distribspec')">Link to the help function</a>
%
%    distribspec generalises the MATLAB function normspec.m for plotting a
%    selected probability density function by shading the portion inside or
%    outside given limits.
%
% Required input arguments:
%
%    pd     : can be one of the following.
%             - A probability density object returned by makedist.
%             - A structure containing two fields:
%               (i) pd.distname, a character array with a valid distribution
%               name (all those accepted by makedist);
%               (ii) pd.x, a numeric vector containing the data sample
%               to be used to estimate the parameters of the distribution.
%             - A structure containing three fields with these features:
%               (i)  pd.distname containing the character array 'user'.
%               (ii) pd.userpdf containing a user defined function handle
%                    expressing the probability density function;
%               (iii) pd.usercdf containing a user defined function handle
%                    expressing the cumulative density function;
%               (iv) pd.x, a numeric vector containing the data sample
%                    to be used to estimate the parameters of the
%                    user-defined distribution. The estimation is done
%                    by maximum likelihood using MATLAB function mle.
%             - A numeric vector containg a sample used to fit a probability
%               distribution object with nonparametric kernel-smoothing.             
%             Other optional fields can be set to help improving the speed
%             and quality of the results obtained for user-defined
%             distributions.
%             In particular, when a data sample is provided in pd.x, the  
%             user can add to structue pd additional fields to control the
%             estimation of the distribution parameters made by distribspec 
%             using matlab function mle:
%             - pd.mleStart, used as initial parameter values. The choice
%               of the starting point is often crucial for the convergence
%               of mle. If the initial parameter values are far from the
%               true maximum likelihood estimate, we can obtain underflow 
%               in the distribution functions and infinite loglikelihoods. 
%               A reasonable choice for the initial value of a location
%               parameter could be pd.mleStart = mean(pd.x). Similarily,
%               for a scale parameter we could set pd.mleStart = std(pd.x).
%             - pd.mleLowerBound, which indicates Lower bounds for
%               distribution parameters. For example, if we trust
%               pd.mleStart, we could set pd.mleLowerBound = pd.mleStart/2.
%             - pd.mleUpperBound, which indicates teh  Upper bounds for
%               distribution parameters. Along the previous case, we could 
%               set pd.mleUpperBound = pd.mleStart*2;
%
%
%    specs  : the lower and upper limits of the shading area. Two element
%             vector. If there is no lower limit, then specs(1)=-Inf. If
%             there is no upper limit, then specs(2)=Inf.
%
%    region : the region to shade. Character array. It can be either the
%             portion 'inside' or 'outside' of the spec limits.
%
%
% Optional input arguments:
%
%    userColor :   The color of the shaded area. Character, 2-element
%                  character array, RGB triplet, 2-row RGB triplet. The
%                  character can be any of the LineSpec colors properties
%                  of MATLAB plot funtion. The RGB triplet specification is
%                  defined as usual, that is with three numbers taking
%                  values in [0 1]: to generate the triplet, the user can
%                  use FSDA function FSColors.
%                   Example - 'userColor', 'r'
%                   Data Types - character
%                   Example - 'userColor', 'rb'
%                   Data Types - 2-character array
%                   Example - 'userColor', FSColors.lightgrey.RGB
%                   Data Types - 3-element array for a RGB triplet
%                   Example - 'userColor', [FSColors.lightgrey.RGB ; FSColors.grey.RGB]
%                   Data Types - 2-row with RGB triplets
%
%   evalPoints :    Evaluation points. Scalar. Number of points where the
%                   the pdf is evaluated and plotted. Default is 1000. The
%                   default is ok for most cases. However, when distribspec
%                   appears too slow, try reducing the number of evaluation
%                   points. If on the contrary the plot looks not
%                   sufficiently smooth, increase the number of evaluation
%                   points.
%                   Example - 'evalPoints', 1
%                   Data Types - scalar
%                   Example - 'evalPoints', [1.5 3.2]
%                   Data Types - array
%
% Output:
%
%    p:   Probability covered by the shaded area. Scalar. It is a value in [0 1].
%
%    h:   Handle to the line objects. Graphic object. Graphic handle.
%
%
% More About: See https://www.mathworks.com/help/stats/fitting-custom-univariate-distributions-part-2.html
%              
%
% See also: normspec, makedist, fitdist, mle
%
% References:
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('distribspec')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Use with makedist, default settings.
    % A standard normal distribution in [-1 1].
    pd     = makedist('Normal');
    specs  = [-1 1];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region);
    pause(0.5);
%}

%{
    %% Use with makedist, region inside.
    % A Gamma with parameter values a = 3 and b = 1, in [2 3].
    close all
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [2 3];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region);
    pause(0.5);
%}

%{
    % Use with makedist, region outside.
    % A Beta with parameter values a = 2 and b = 4, outside [0.2 0.4].
    close all
    pd = makedist('Beta','a',2,'b',4);
    specs  = [0.2 0.4];
    region = 'outside';
    [p, h] = distribspec(pd, specs, region);
    pause(0.5);
%}

%{
    % Use with makedist, specs with inf.
    % A Beta with parameter values a = 2 and b = 4, in [0.4 inf].
    close all
    pd = makedist('Beta','a',2,'b',4);
    specs  = [0.4 inf];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region);
    pause(0.5);
%}

%{
    % Use with makedist, specs with -inf.
    
    % A Gamma with parameter values a = 3 and b = 1, up to 2.
    close all
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [-inf 2];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region);

    % Without specification of region the default is inside, therefore,
    % the output is identical
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [-inf 2];
    [p, h] = distribspec(pd, specs);
    
    cascade;
    pause(0.5);
%}

%{
    % Use with makedist, no specs.
    % A Gamma as above, without specification of specs: returns an error.
    close all
    pd = makedist('Gamma','a',3,'b',1);
    region = 'inside';
    try
        [p, h] = distribspec(pd, [], region);
    catch
        disp(sprintf('%s\n\n%s\n%s',...
            'Without specification of specs, distribspec returned this error:' , ...
            'Error using distribspec' , ...
            'Please provide valid lower and upper limits for the area to highlight'));
    end
    pause(0.5);
%}

%{
    % Use with makedist, specs outside support.
    % A Gamma with parameter values a = 3 and b = 1, up to -1. 
    % Nothing is colored.
    close all
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [-inf -1];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region);
    pause(0.5);
%}

%{
    % Use with makedist, using userColor.

    % A Gamma as above, using userColor with standard one-character
    % specification.
    close all
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [-inf 2];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region, 'userColor','r');
    
    % same, but using userColor with RGB triplet specification.
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [-inf 2];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region, 'userColor',[1 0.5 0.5]);
    
    cascade;
    pause(0.5);
%}

%{
    %% Use with makedist, using userColor for each outside region.

    close all
    pd = makedist('HalfNormal','mu',0,'sigma',1.5);
    specs  = [1 3];
    region = 'outside';
    useColor = 'cg';
    [p, h] = distribspec(pd, specs, region, 'userColor', useColor);

    useColor = [1 0.5 0.5 ; 0.5 0.8 0.3];
    [p, h] = distribspec(pd, specs, region, 'userColor', useColor);

    cascade;
    pause(0.5);
%}

%{
    % Use with makedist, using userColor with RGB triplet. 

    % A Half Normal using userColor with RGB triplet specification.
    % The triplet is returned by FSColors.
    close all
    pd = makedist('HalfNormal','mu',0,'sigma',1.5);
    specs  = [-inf 1];
    region = 'inside';
    RGB_vector = FSColors.lightgrey.RGB;
    [p, h] = distribspec(pd, specs, region, 'userColor',RGB_vector);
    pause(0.5);
%}

%{
    %% Use with distname. 

    % A sample of n=100 elements extracted from a Nakagami.
    close all
    rng(12345);
    distname    = 'Nakagami';
    x           = random(distname,0.5,1,[100,1]);
    pd          = struct;
    pd.x        = x;
    pd.distname = distname;
    specs  = [-inf 1];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region);
    pause(0.5);
%}

%{
    %% Use with a sample only. 

    % A sample of n=100 elements extracted from a T(5) without distname.
    close all
    rng(12345);
    x      = random('T',5,[100,1]);
    specs  = [-inf 1.5];
    region = 'inside';
    [p, h] = distribspec(x, specs, region);
    pause(0.5);
%}

%{ 
    % Sample from non-parametric used-defined distribution.

    close all
    % Data from user-defined distribution, Benford (non-parametric)
    % no parameters to estimate

    userpdf = @(n) log10(1 + 1./n);
    usercdf = @(n) log10(1 + n);

    % random sample from Benford
    n = 5000;
    u = rand(n,1);
    X = 10.^u;
    
    % Apply distribspec with userpdf
    specs  = [-inf 5];
    region = 'inside';
    pd          = struct;
    pd.x        = X;
    pd.distname = 'user';
    pd.userpdf  = userpdf;

    [p, h] = distribspec(pd, specs, region);   
    pause(0.5);
%}

%{
    %% Sample from parametric used-defined distribution.

    % Data from user-defined distribution, negative exponential with one parameter 
    close all
    lambda = 5;
    userpdf = @(data,lambda) lambda*exp(-lambda*data);
    usercdf = @(data,lambda) 1-exp(-lambda*data);

    rng(12345);
    n = 500;
    X = zeros(n,1);
    u = rand(n,1);
    for i = 1:numel(u)
      fun  = @(x)integral(@(x)userpdf(x,lambda),eps,x) - u(i);
      X(i) = fzero(fun,0.5);
      % The previous two lines have the solution below, but exemplify the
      % approach for a generic function without closed formula.
      % X(i) = -(1/lambda)*log(1-u(i)); % p. 211 Mood Graybill and Boes
    end
    
    %Estimate the parameter lambda of the custom distribution.
    parhat = mle(X,'pdf',userpdf,'cdf',usercdf,'start',0.05);
    
    % plot data as a normalised histogram and the mle of the pdf
    histogram(X,100,'Normalization','pdf');
    hold on
    plot(X,userpdf(X,parhat),'.')
    title({['The sample generated by ' char(userpdf)] , ['with $\lambda=$' num2str(lambda)]} , 'Interpreter' , 'Latex');
    hold off

    % Use distribspec with both usercdf and userpdf
    specs  = [-inf 0.2];
    region = 'inside';
    pd          = struct;
    pd.x        = X;
    pd.distname = 'user';
    pd.usercdf  = usercdf;
    pd.userpdf  = userpdf;
    [p, h] = distribspec(pd, specs, region);

    % now with userpdf
    pd          = struct;
    pd.x        = X;
    pd.distname = 'user';
    pd.userpdf  = userpdf;
    [p, h] = distribspec(pd, specs, region);

    % and now with usercdf
    pd          = struct;
    pd.x        = X;
    pd.distname = 'user';
    pd.usercdf  = usercdf;
    [p, h] = distribspec(pd, specs, region);

    cascade;
    pause(0.5);
%}

%{
    %% User-defined function, with reduced number of evaluation points.

    % Data from user-defined distribution: Pareto (two parameters). This
    % takes a while to complete, because of difficult integration (the
    % function has singularities). Therefore, we reduce the ealuation
    % points only to 50. In this case, being the function very smooth in
    % the region of interest, the quality is not affected.
    
    close all
    alpha0  = 2; 
    xm0     = 4;
    
    % data
    rng(12345);
    n = 500;
    X = zeros(n,1);
    u = rand(n,1);
    for i = 1:numel(u)
      X(i) = alpha0*(1-u(i)).^(-1/xm0);
    end

    % heaviside function
    hvsd = @(y) [0.5*(y == 0) + (y > 0)];
    userpdf0 = @(x) hvsd(x-xm0) .* (alpha0*4^alpha0) ./ (x.^(alpha0+1));
    fplot(userpdf0);
    xlim([-5 20]);
    title({'Pareto distribution' , '$ hvsd(x-xm_{0}) \cdot (\alpha_{0} \cdot 4^{\alpha_{0}}) / (x^{\alpha_{0}+1})$' } , 'Interpreter','latex' , 'Fontsize' , 20);
    pause(0.5);

    userpdf = @(x, alpha) hvsd(x-4) .* ((alpha*4^alpha) ./ (x.^(alpha+1)));
    %userpdf = @(x, alpha, xm) hvsd(x-xm) .* ((alpha*xm^alpha) ./ (x.^(alpha+1)));
    %usercdf = @(x, alpha, xm) hvsd(x-xm) .* (1 - (xm./x).^alpha);

    % Apply distribspec with userpdf
    specs  = [5 10];
    region = 'inside';
    pd          = struct;
    pd.x        = X;
    pd.distname = 'user';
    pd.userpdf  = userpdf;
    pd.mleStart = mean(pd.x);
    pd.mleLowerBound = pd.mleStart/2;
    pd.mleUpperBound = pd.mleStart*2;

    [p, h] = distribspec(pd, specs, region, 'evalPoints', 50);

    cascade;
    pause(0.5);
%}

%{
    %% User-defined function starting from a default matlab distribution.

    % Here the user takes one of the functions in makedist, possibly
    % changing the number of parameters.
    
    close all
    alpha0  = 2; 
    xm0     = 4;

    X = wblrnd(1,1,[1000,1]) + 10;
    histogram(X,'Normalization','pdf');
    title({'A sample from the weibull' , '$a=1, b=1, c=10$'},'Fontsize',20,'Interpreter','latex');

    % A Weibull distribution with an extra parameter c
    userpdf = @(x,a,b,c) wblpdf(x-c,a,b);

    % Apply distribspec with userpdf
    specs  = [11 12];
    region = 'inside';
    pd          = struct;
    pd.x        = X;
    pd.distname = 'user';
    pd.userpdf  = userpdf;
    pd.mleStart = [5 5 5];
    pd.mleLowerBound = [0 0 -Inf];
    pd.mleUpperBound = [Inf Inf min(X)];

    % MLE parameters. The scale and shape parameters must be positive,
    % and the location parameter must be smaller than the minimum of the
    % sample data.
    
    [p, h] = distribspec(pd, specs, region, 'evalPoints', 100);
    cascade
    pause(0.5);
%}


%% Beginning of code

fittedUsingKernel   = false;
fittedUsingDistname = false;
userpdf             = [];
usercdf             = [];

%% Main parameters

if nargin<3 || isempty(region)
    region='inside';
end

if numel(specs) ~= 2 || ~isnumeric(specs)
    error('FSDA:distribspec:BadSpecsSize','Please provide valid lower and upper limits for the area to highlight');
end

if ~strcmp(region,'inside') && ~strcmp(region,'outside')
    error('FSDA:distribspec:BadRegion','The region can be either "inside" or "outside"');
end

%% Optional parameters

options = struct('userColor', 'b', 'evalPoints', 1000);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)

    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end

    % Check if the specified optional arguments were present in structure
    % options
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end

end
% Write in structure 'options' the options chosen by the user
for i=1:2:length(varargin)
    options.(varargin{i})=varargin{i+1};
end

% Assign the values for the optional arguments
userColor  = options.userColor;
evalPoints = options.evalPoints;

% Evaluation points
prob = linspace(0.001,0.999,evalPoints)';

%% Checks on the structure provided by the user

if isstruct(pd)

    if ~and(isfield(pd,'distname') , isfield(pd,'x'))
        error('FSDA:distribspec:BadFilds','Bad filds or filed names: please specify "distname" and "x".');
    end
    if ~any( strcmp( strtrim(makedist),strtrim(lower(pd.distname)) ) )
        if ~strcmp('user' , strtrim(lower(pd.distname)) )
            error('FSDA:distribspec:BadDistribFildNames','Bad distribution name: see makedist for a comprehensive list.');
        else
            if ~isfield(pd,'userpdf') && ~isfield(pd,'usercdf')
                error('FSDA:distribspec:BadUserDistribFild','User distribution missing: please specify your distribution as function handle.');
            else
                % check if there is a user-specific pdf and/or cdf function
                if isfield(pd,'userpdf')
                    userpdf   = pd.userpdf;
                    parnum_p  = nargin(userpdf)-1;
                else
                    parnum_p  = -1; % there is no userpdf
                end
                if isfield(pd,'usercdf')
                    usercdf = pd.usercdf;
                    parnum_c  = nargin(usercdf)-1;
                else
                    parnum_c  = -1; % there is no usercdf
                end
                % parameters of the cdf/pdf
                parnum = max(parnum_c,parnum_p) ;
                
                if ~isfield(pd,'mleStart')
                    pd.mleStart = repmat(0.5,parnum,1);
                end
                if ~isfield(pd,'mleLowerBound')
                    pd.mleLowerBound = [];
                end
                if ~isfield(pd,'mleUpperBound')
                    pd.mleUpperBound = [];
                end
            end
        end
    end
    if ~isnumeric(pd.x)
        error('FSDA:distribspec:BadSampleX','Bad sample: the array x must be numeric.');
    end
end

isuserdistrib = false;

%% Depending on pd, computation of x and y

% In the next 'if' statement:
% - we check inputs and invalid arguments
% - compute the x and y values for predefined quantiles and the
%   selected distribution/data/parameters

if isobject(pd) && isprop(pd,'DistributionName')

    % CASE 1. pd CONTAINS A PROBABILITY DISTRIBUTION OF makedist.m FUNCTION

    % set of values based on the given probability density object created
    % by the user outside distribspec
    x = icdf(pd,prob);
    y = pdf(pd, x);

    % plot(x,y,'-')

elseif isstruct(pd)

    % CASE 2. pd CONTAINS A PROBABILITY DISTRIBUTION AND A SAMPLE

    if isempty(userpdf) && isempty(usercdf)

        % CASE 2.1: THE PROBABILITY DISTRIBUTION IS ONE OF makedist

        fittedUsingDistname = true;

        % set the values based on the given probability density name and sample
        sample = pd.x; sample = sample(:); % fitdist requires x to be a column vector
        pd = fitdist(sample,pd.distname);
        x = icdf(pd,prob);
        y = pdf(pd, x);
        
        % plot(x,y,'-')

    else

        % CASE 2.2: THE PROBABILITY DISTRIBUTION IS USER-DEFINED
        
        isuserdistrib = true;
        mleOpt        = statset('FunValCheck','off');
        %fminuncOpt    = optimset('Display','off','MaxIter',10000,'TolX',10^-30,'TolFun',10^-30);

        % The available sample
        sample = pd.x;
        sample = sample(:);

        % Parameter estimation
        if parnum_c > 0 && parnum_p > 0
            % Both pdf and cdf are specified
            parhat = mle(sample,'pdf',userpdf,'cdf',usercdf,'start',pd.mleStart,'LowerBound',pd.mleLowerBound,'UpperBound',pd.mleUpperBound,'Options',mleOpt);

        elseif parnum_p > 0 && parnum_c == -1
            % The pdf is specified
            parhat   = mle(sample,'pdf',userpdf,'start',pd.mleStart,'LowerBound',pd.mleLowerBound,'UpperBound',pd.mleUpperBound,'Options',mleOpt);

        elseif parnum_c > 0 && parnum_p == -1
            % The cdf is specified (mle does not work)

            % workaround to find the icdf
            if min(sample)>=0
                xx = ksdensity(sample,sample,'Function','icdf','Support','Positive');
            else
                xx = ksdensity(sample,sample,'Function','icdf');
            end
            ii = xx>0;

            % now find the parameter(s) by minimizing the rmse
            rmse   = @(xpar)sqrt(sum((sample(ii) - usercdf(xx(ii),xpar)).^ 2));
            parhat = fminsearch(rmse,0);
            %options = optimset('Display','iter','PlotFcns',@optimplotfval);
            %parhat = fminsearch(rmse,0,options);

        end

        % Computation of the x and y values 
        if parnum_c > 0 && parnum_p > 0
            % The user has provided both the pdf and cdf

            % Compute the ICDF from the CDF
            switch parnum
                case 1
                    userICDF = @(p) fzero( @(pp)usercdf(pp,parhat) - p , 0);
                case 2
                    userICDF = @(p) fzero( @(pp)usercdf(pp,parhat(1),parhat(2)) - p , 0);
                case 3
                    userICDF = @(p) fzero( @(pp)usercdf(pp,parhat(1),parhat(2),parhat(3)) - p , 0);
            end
            x = arrayfun(userICDF,prob);
            switch parnum
                case 1
                    y = userpdf(x,parhat);
                case 2
                    y = userpdf(x,parhat(1),parhat(2));
                case 3
                    y = userpdf(x,parhat(1),parhat(2),parhat(3));
            end

            % arrayfun replaces the following loop:
            %x = zeros(numel(prob),1);
            %for i=1:numel(prob)
            %    x(i) = userICDF(prob(i),parhat);
            %end

            % plot(x,y,'-')

        elseif parnum_p > 0 && parnum_c == -1
            % the user has provided the pdf

            % Compute the CDF from the userpdf
            switch parnum
                case 1
                    userCDF = @(x)integral(@(x)userpdf(x,parhat),eps,x,'ArrayValued',1);
                case 2
                    userCDF = @(x)integral(@(x)userpdf(x,parhat(1),parhat(2)),eps,x,'ArrayValued',1);
                case 3
                    userCDF = @(x)integral(@(x)userpdf(x,parhat(1),parhat(2),parhat(3)),eps,x,'ArrayValued',1);
            end
            % Compute the ICDF from the CDF 
            initialGuess = eps; % check first if initialGuess leads to a finite value for the integrand 
            fzeroOpt = optimset('MaxFunEvals',100,'MaxIter',100);
            userICDF = @(p) fzero( @(pp)userCDF(pp) - p , initialGuess, fzeroOpt);
            x = arrayfun(userICDF,prob);
            switch parnum
                case 1
                    y = userpdf(x,parhat);
                case 2
                    y = userpdf(x,parhat(1),parhat(2));
                case 3
                    y = userpdf(x,parhat(1),parhat(2),parhat(3));
            end

            % plot(x,y,'-')

        elseif parnum_c > 0 && parnum_p == -1
            % the user has provided the cdf

            initialGuess = 0; % check first if initialGuess leads to a finite value for the integrand
            switch parnum
                case 1
                    userICDF = @(p) fzero( @(pp)usercdf(pp,parhat) - p , initialGuess);
                case 2
                    userICDF = @(p) fzero( @(pp)usercdf(pp,parhat(1),parhat(2)) - p , initialGuess);
                case 3
                    userICDF = @(p) fzero( @(pp)usercdf(pp,parhat(1),parhat(2),parhat(3)) - p , initialGuess);
            end
                    
            x = arrayfun(userICDF,prob);

            % Workaround for the above
            % if min(sample)>=0
            %     xx = ksdensity(sample,linspace(min(sample),max(sample),numel(prob)),'Function','icdf','Support','Positive');
            % else
            %     xx = ksdensity(sample,linspace(min(sample),max(sample),numel(prob)),'Function','icdf');
            % end
            % ii = xx>0; x  = xx(ii); x  = x(:);

            % Compute the pdf values as d(CDF) / d(x)
            h = 1e-9;
            %complex step differention: imag(f(x+h*1i)/h);
            switch parnum
                case 1
                    y = imag(usercdf(prob+h*1i,parhat) / h);
                case 2
                    y = imag(usercdf(prob+h*1i,parhat(1),parhat(2)) / h);
                case 3
                    y = imag(usercdf(prob+h*1i,parhat(1),parhat(2),parhat(3)) / h);
            end

            % there are other ways to differentiate:
            %y0 = diff([eps ; usercdf(prob,parhat)]) ./ diff([eps ; prob]);
            %y1 = gradient(usercdf(prob,parhat)) ./ gradient(prob);
            %y2 = gradient(usercdf(prob(ii),parhat),prob(ii));
            % finite difference (f(x) - f(x-h)) / h;
            %y3 = (usercdf(prob(ii),parhat) - usercdf(prob(ii)-h,parhat)) / h
            % central difference (f(x+h)-f(x-h))/(2*h);
            %y4 = (usercdf(prob(ii)+h,parhat) - usercdf(prob(ii)-h,parhat)) / (2*h);

            %[yy,xx] = hist(x); plot(xx,yy/sum(pp)); %PDF from sample
            %[ff,xx] = ecdf(x); plot(xx,ff);         %CDF from sample

            % plot(x,y,'-') 

        else
            % the user distribution has no parameters to estimate

            % Workaround for the above
            if min(sample)>=0
                x = ksdensity(sample,'Function','icdf','NumPoints',numel(prob),'Support','Positive');
            else
                x = ksdensity(sample,'Function','icdf','NumPoints',numel(prob));
            end
            ii = x>0; 
            x = x(ii);
            x = x(:);

            y = userpdf(x);

            % % Compute the CDF from the userpdf: we use the compact version below
            % usercdf = @(x)integral(@(x)userpdf(x),eps,x,'ArrayValued',1);
            % 
            % % Compute the ICDF from the CDF, as in the previous case
            % userICDF = @(p) fzero( @(pp)usercdf(pp) - p , 0.5);
            % 
            % x = zeros(numel(prob),1);
            % for i=1:numel(prob)
            %     x(i) = userICDF(prob(i));
            % end
            % y = userpdf(x);

            % plot(x,y,'-')
        end

        % normalization of the pdf
        %area = trapz(y);
        %y = y/area;

        % workaround to obtain a matlab distribution object, used for the plot
        if min(x)>=0 
            pd = fitdist(x,'Kernel','Support','Positive');
        else
            pd = fitdist(x,'Kernel');
        end

    end

elseif isnumeric(pd)

    %CASE 3. pd CONTAINS JUST A DATA SAMPLE

    fittedUsingKernel = true;

    x = pd;
    if min(x)>=0
        pd = fitdist(x,'Kernel','Support','Positive');
    else
        pd = fitdist(x,'Kernel');
    end

    x = icdf(pd,prob);
    y = pdf(pd, x);

else
    error('FSDA:distribspec:BadDistribOption','"pd" option has been mispecified: read carefully the help.');
end


%% set upper and lower bound

if nargin<2 || isempty(specs)
    specs=[min(x) max(x)];
    emptyspecs = true;
else
    emptyspecs = false;
end

% swap the specs if they are reversed
if specs(1) > specs(2)
    specs = fliplr(specs);
end
lb    = specs(1);
ub    = specs(2);
lbinf = isinf(lb);
ubinf = isinf(ub);

if lbinf && ubinf
    error('FSDA:distribspec:BadSpecsInfinite','Both limits are infinite');
end


%% compute p

if strcmp(region,'outside')
    if emptyspecs
        p=0;
        strprob = ['No specs given by the user: probability is ' num2str(p)];
    else
        if lbinf
            p = 1-cdf(pd,ub); % P(t > ub)
            strprob = ['Probability greater than upper bound is ' num2str(p)];
        elseif ubinf
            p = cdf(pd,lb); % P(t < lb)
            strprob = ['Probability smaller than lower bound is ' num2str(p)];
        else
            p = cdf(pd,lb) + (1-cdf(pd,ub)); % P(t < lb) + Pr(t > ub)
            strprob = ['The outside region $P(t < lb) + Pr(t > ub)$ is ' num2str(p)];
        end
    end
else  % if strcmp(region,'inside')
    if emptyspecs
        p=1;
        strprob = ['No specs given by the user: probability is ' num2str(p)];
    else
        if lbinf
            p = cdf(pd,ub); % P(t < ub)
            strprob = ['Probability lower than upper bound is ' num2str(p)];
        elseif ubinf
            p = 1-cdf(pd,lb); % P(t > lb)
            strprob = ['Probability greater than lower bound is ' num2str(p)];
        else
            p = cdf(pd,ub) - cdf(pd,lb);  % P(lb < t < ub)
            strprob = ['The inside region $P(lb < t < ub)$ is ' num2str(p)];
        end
    end
end


%% plot the distribution

% initialise figure
nspecfig = figure;
nspecaxes = axes;
set(nspecaxes, 'Parent', nspecfig);
set(nspecaxes,'Nextplot','add');
hh    = plot(x,y,'k-','LineWidth',1.5);

if exist('sample','var')
    xlim([min(x) , (max(x)+max(sample))/2]);
end
xlims = get(nspecaxes,'Xlim');

% compute the endpoints of the spec limit lines and plot limit lines
pll =  [max(xlims(1),min(x));max(xlims(1),min(x))];
ypll = [0;eps];
if lbinf
    ll =  pll;
    yll = ypll;
else
    ll =  [lb; lb];
    if isuserdistrib
        switch parnum
            case 0
                yll = [0; userpdf(lb)];
            case 1
                yll = [0; userpdf(lb,parhat)]; 
            case 2
                yll = [0; userpdf(lb,parhat(1),parhat(2))]; 
            case 3
                yll = [0; userpdf(lb,parhat(1),parhat(2),parhat(3))]; 
        end
    else
        yll = [0; pdf(pd, lb)];  
    end
end

pul =  [xlims(2);xlims(2)];
ypul = [eps;0];
if ubinf
    ul =  pul;
    yul = ypul;
else
    ul  = [ub; ub];
    if isuserdistrib
        if parnum==0
            yul = [userpdf(ub); 0];  
        else
            if parnum_p > 0
                % the user has provided a pdf
                switch parnum
                    case 1
                        yul = [userpdf(ub,parhat); 0]; 
                    case 2
                        yul = [userpdf(ub,parhat(1),parhat(2)); 0];                         
                    case 3
                        yul = [userpdf(ub,parhat(1),parhat(2),parhat(3)); 0]; 
                end
            else
                % the user has provided a cdf only
                [~, indexOfMinDist] = min(abs(x-ub));
                yul = [y(indexOfMinDist); 0]; 
            end
        end
    else
        yul = [pdf(pd, ub); 0];  
    end
end

% fill the shade area
switch region
    case 'inside'
        if ubinf
            k = find(x > lb);
            hh1 = plot(ll,yll,'b-');
        elseif lbinf
            k = find(x < ub);
            hh1 = plot(ul,yul,'b-');
        else
            k = find(x > lb & x < ub);
            hh1 = plot(ll,yll,'b-',ul,yul,'b-');
        end
        xfill = [ll; x(k); ul];
        yfill = [yll; y(k); yul];

        if (ischar(userColor) && numel(userColor)>1)
            userColor = userColor(1);
            warning('FSDA:distribspec:BaduserColor','userColor has more than one color specification: only the first one is used for the inside region');
        end
        if (isnumeric(userColor) && size(userColor,1)>1)
            userColor = userColor(1,:);
            warning('FSDA:distribspec:BaduserColor','userColor has more than one color specification: only the first one is used for the inside region');
        end

        fill(xfill,yfill,userColor);

    case 'outside'
        if ubinf
            k1 = find(x < lb);
            k2=[];
            hh1 = plot(ll,yll,'b-');
        elseif lbinf
            k1=[];
            k2 = find(x > ub);
            hh1 = plot(ul,yul,'b-');
        else
            k1 = find(x < lb );
            k2=find(x > ub);
            hh1 = plot(ll,yll,'b-',ul,yul,'b-');
        end

        % fill regions with user-defined or default color
        if (ischar(userColor) && numel(userColor)==1) || (isnumeric(userColor) && size(userColor,1)==1)
            xfill = [pll;  x(k1); ll          ; ul;          x(k2); pul  ];
            yfill = [ypll; y(k1); flipud(yll) ; flipud(yul); y(k2); ypul ];
            fill(xfill,yfill,userColor);
        elseif (ischar(userColor) && numel(userColor)==2) || (isnumeric(userColor) && size(userColor,1)==2)
            xfill1 = [pll;  x(k1); ll          ];
            yfill1 = [ypll; y(k1); flipud(yll) ];
            xfill2 = [ul;          x(k2); pul  ];
            yfill2 = [flipud(yul); y(k2); ypul ];
            if ischar(userColor)
                fill(xfill1,yfill1,userColor(1));
                fill(xfill2,yfill2,userColor(2));
            else
                fill(xfill1,yfill1,userColor(1,:));
                fill(xfill2,yfill2,userColor(2,:));
            end
        else
            warning('FSDA:distribspec:BaduserColor','userColor is wrong: a default is used');
            xfill = [pll;  x(k1); ll          ; ul;          x(k2); pul  ];
            yfill = [ypll; y(k1); flipud(yll) ; flipud(yul); y(k2); ypul ];
            fill(xfill,yfill,'b');
        end
end

%%  Add title and labels to the plot and return the handles

if fittedUsingKernel
    % CASE 3: there is no distribution specified
    title({'Nonparametric kernel-smoothing distribution' , 'Fit by fitdist.m with default options'}, 'interpreter' , 'latex' , 'FontSize', 14);

elseif ~isempty(userpdf)
    % CASE 2.2.a: this is for a user-defined distribution: pdf
    if parnum>0
        strpar = ['with mle parameter ' , '$\hat{\theta}$' '=' num2str(parhat) ' '];
    else
        strpar = '';
    end
    title({'User-defined PDF' ; ['$' char(userpdf) '$'] ; strpar ; strprob ; ['$lb =$ ' num2str(lb) ' -- ' '$ub =$ ' num2str(ub)]}, 'interpreter' , 'latex' , 'FontSize', 14);

elseif ~isempty(usercdf)
    % CASE 2.2.b: this is for a user-defined distribution: cdf
    if parnum>0
        strpar = ['with argmin(rmse) parameter ' , '$\hat{\theta}$' '=' num2str(parhat) ' '];
    else
        strpar = '';
    end
    title({'User-defined CDF' ; ['$' char(usercdf) '$'] ; strpar ; strprob ; ['$lb =$ ' num2str(lb) ' -- ' '$ub =$ ' num2str(ub)]}, 'interpreter' , 'latex' , 'FontSize', 14);

else
    % CASE 1 & 2.1: this is if the user has selected one of the matlab distributions
    numpar = length(pd.ParameterNames);
    strpar = [pd.DistributionName ' with '];
    for np=1:numpar
        if fittedUsingDistname
            strpar = [strpar , '$\hat{' char(pd.ParameterNames(np)) '}$' '=' num2str(pd.ParameterValues(np)) ' ']; %#ok<AGROW>
        else
            strpar = [strpar , char(pd.ParameterNames(np)) '=' num2str(pd.ParameterValues(np)) ' ']; %#ok<AGROW>
        end
    end
    title({strpar ; strprob ; ['$lb =$ ' num2str(lb) ' -- ' '$ub =$ ' num2str(ub)]}, 'interpreter' , 'latex' , 'FontSize', 14);
end

% axis and labels
xaxis = refline(0,0);
set(xaxis,'Color','k');
ylabel('Density value', 'interpreter' , 'latex' , 'FontSize', 12);
xlabel('Critical value', 'interpreter' , 'latex' , 'FontSize', 12);

if nargout > 1
    h = [hh; hh1];
end

end
%FScategory:UTISTAT
