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
% Output:
%
%    p:   Probability covered by the shaded area. Scalar. It is a value in [0 1].
%
%    h:   Handle to the line objects. Graphic object.
%
%
% Optional Output:
%
%
% References:
%
%
% See also: normspec, makedist, fitdist
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
    % A standard normal distribution in [-1 1].
    pd     = makedist('Normal','mu',0,'sigma',1);
    specs  = [-1 1];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region);
%}

%{
    %% A Gamma with parameter values a = 3 and b = 1, in [2 3].
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [2 3];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region);
%}

%{
    %% A Gamma with parameter values a = 3 and b = 1, outside [2 3].
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [2 3];
    region = 'outside';
    [p, h] = distribspec(pd, specs, region);
%}

%{
    % A Gamma with parameter values a = 3 and b = 1, in [2 inf].
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [2 inf];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region);
%}

%{
    % A Gamma with parameter values a = 3 and b = 1, up to 2.
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [-inf 2];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region);
%}

%{
    % A Gamma as above, without specification of region (default is
    % inside).
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [-inf 2];
    [p, h] = distribspec(pd, specs);
%}

%{
    % A Gamma as above, without specification of specs.
    pd = makedist('Gamma','a',3,'b',1);
    region = 'inside';
    [p, h] = distribspec(pd, [], region);
%}

%{
    % A Gamma as above, without specification of specs.
    pd = makedist('Gamma','a',3,'b',1);
    region = 'outside';
    [p, h] = distribspec(pd, [], region);
%}

%{
    % A Gamma with parameter values a = 3 and b = 1, up to -1.
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [-inf -1];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region);
%}

%{
    % A Gamma as above, using userColor with standard one-character
    % specification.
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [-inf 2];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region, 'userColor','r');
%}

%{
    % A Gamma as above, using userColor with RGB triplet specification.
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [-inf 2];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region, 'userColor',[1 0.5 0.5]);
%}

%{
    %% A Gamma using userColor with a color for each outside region
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [2 3];
    region = 'outside';
    useColor = 'cg';
    [p, h] = distribspec(pd, specs, region, 'userColor', useColor);

    useColor = [1 0.5 0.5 ; 0.5 0.8 0.3];
    [p, h] = distribspec(pd, specs, region, 'userColor', useColor);

    cascade;
%}

%{
    %% A Gamma as above, using userColor with RGB triplet specification
    % returned by FSColors.
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [-inf 2];
    region = 'inside';
    RGB_vector = FSColors.lightgrey.RGB;
    [p, h] = distribspec(pd, specs, region, 'userColor',RGB_vector);
%}

%{
    %% A sample of n=100 elements extracted from a Gamma, with distname.
    rng(12345);
    distname    = 'Gamma';
    x           = random(distname,3,1,[100,1]);
    pd          = struct;
    pd.x        = x;
    pd.distname = distname;
    specs  = [-inf 2];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region);
%}

%{
    %% A sample of n=100 elements extracted from a T(5) without distname.
    rng(12345);
    x      = random('T',5,[100,1]);
    specs  = [-inf 2];
    region = 'inside';
    [p, h] = distribspec(x, specs, region);
%}

%{
    %% A sample from a user-defined distribution with given parameter lambda. 
    lambda = 5;
    userpdf = @(data,lambda) lambda*exp(-lambda*data);
    usercdf = @(data,lambda) 1-exp(-lambda*data);

    rng(12345);
    n = 500;
    X = zeros(n,1);
    u = rand(n,1);
    for i = 1:numel(u)
      fun  = @(x)integral(@(x)userpdf(x,lambda),0,x) - u(i);
      X(i) = fzero(fun,0.5);
      % The previous two lines have the solution below, but exemplify the
      % approach for a generic function without closed formula.
      % X(i) = -(1/lambda)*log(1-u(i)); % p. 211 Mood Graybill and Boes
    end
    
    %Estimate the parameter, lambda, of the custom distribution for the censored sample data.
    parhat = mle(X,'pdf',userpdf,'cdf',usercdf,'start',0.05);
    
    % plot data as a normalised histogram and the mle of the pdf
    histogram(X,100,'Normalization','pdf');
    hold on
    plot(X,userpdf(X,parhat),'.')
    title({['The sample generated by ' char(userpdf)] , ['with $\lambda=$' num2str(lambda)]} , 'Interpreter' , 'Latex');

    % now, use distribspec with userpdf
    specs  = [-inf 0.1];
    region = 'inside';
    pd          = struct;
    pd.x        = X;
    pd.distname = 'user';
    pd.usercdf  = usercdf;

    [p, h] = distribspec(pd, specs, region);

    rmfield(pd,'usercdf');
    pd.userpdf  = userpdf;
    [p, h] = distribspec(pd, specs, region);

    cascade;
%}

%% Beginning of code

fittedUsingKernel   = false;
fittedUsingDistname = false; 
userpdf             = [];
usercdf             = [];

% In the next 'if' statement: 
% - we check inputs and invalid arguments
% - we generate a set of quantiles and the respective x and y values for
%   the selected distribution/data/parameters

if isobject(pd) && isprop(pd,'DistributionName')
    % CASE 1. pd CONTAINS A PROBABILITY DISTRIBUTION OF makedist.m FUNCTION

    % set of values based on the given probability density object
    prob = (0.0001:0.0004:0.9999)';
    x = icdf(pd,prob);
    y = pdf(pd, x);

elseif isstruct(pd)
    % CASE 2. pd CONTAINS A PROBABILITY DISTRIBUTION AND A SAMPLE
    
    % do the necessary checks on the structure provided by the user
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
            end
        end
    end
    if ~isnumeric(pd.x)
        error('FSDA:distribspec:BadSampleX','Bad sample: the array x must be numeric.');
    end

    if isempty(userpdf) && isempty(usercdf)
        % CASE 2.1: THE PROBABILITY DISTRIBUTION IS ONE OF makedist
        
        fittedUsingDistname = true; 
    
        % set the values based on the given probability density name and sample
        sample = pd.x; sample = sample(:); % fitdist requires x to be a column vector
        pd = fitdist(sample,pd.distname);
        prob = (0.0001:0.0004:0.9999)';
        x = icdf(pd,prob);
        y = pdf(pd, x);

    else
        % CASE 2.2: THE PROBABILITY DISTRIBUTION IS USER-DEFINED 

        sample  = pd.x; sample = sample(:);
        prob    = (0.0001:0.0004:0.9999)';
        parnum  = max(parnum_c,parnum_p) ; %nargin(userpdf)-1; 

        if parnum_c > 0 && parnum_p > 0
            % the user has provided both the pdf and cdf
            parhat = mle(sample,'pdf',userpdf,'cdf',usercdf,'start',0.05);
            x = usercdf(prob,parhat);
            y = userpdf(x,parhat);

        elseif parnum_p > 0 && parnum_c == -1
            % the user has provided the pdf
            parhat  = mle(sample,'pdf',userpdf,'start',0.05);
            usercdf = @(x)integral(@(x)userpdf(x,parhat),0,x,'ArrayValued',1);
            x = zeros(numel(prob),1);
            for i=1:numel(prob)
                x(i) = usercdf(prob(i));
            end
            y = userpdf(x,parhat);

        elseif parnum_c > 0 && parnum_p == -1
            % the user has provided the cdf

            % the following statement does not work: uses the normal
            % parhat = mle(sample,'cdf',usercdf,'start',0.05);
            
            % workaround to find the values to be given to the cdf
            pd = fitdist(sample,'Kernel','Support','positive');
            xx = icdf(pd,sample); 
            ii = find(xx>0); % some values seem badly estimated
            % now find the parameter(s) by minimizing the rmse
            rmse = @(xpar)sqrt(sum((sample(ii) - usercdf(xx(ii),xpar)).^ 2));
            parhat = fminsearch(rmse,0);

            %{  
                % this might be better of the above with some tunings, but
                % it requires the optimization toolbox 
                parhat = lsqcurvefit(usercdf,0,xx(ii),sample(ii));
            %}

            x = zeros(numel(prob),1);
            for i=1:numel(prob)
                x(i) = usercdf(prob(i),parhat);
            end

            y = diff(x) ./ diff(prob);
            x = prob(1:end-1) + diff(prob)./2;

            %  plotyy(prob,usercdf(prob,5),x,y);

        else
            % the user distribution has no parameters to estimate
            x = usercdf(prob);
            y = userpdf(x);
        end

        % workaround to obtain a matlab distribution object, used for the plot
        pd = fitdist(x,'Kernel','Support','positive');
        
    end

elseif isnumeric(pd)
    %CASE 3. pd CONTAINS JUST A DATA SAMPLE

    fittedUsingKernel = true;

    x = pd;
    pd = fitdist(x,'Kernel','Support','positive');
    prob = (0.0001:0.0004:0.9999)';
    x = icdf(pd,prob);
    y = pdf(pd, x);

else
    error('FSDA:distribspec:BadDistribOption','"pd" option has been mispecified: read carefully the help.');
end

if nargin<3 || isempty(region)
    region='inside';
end

if nargin<2 || isempty(specs)
    specs=[min(x) max(x)];
    emptyspecs = true;
else
    emptyspecs = false;
end

if numel(specs) ~= 2 || ~isnumeric(specs)
    error('FSDA:distribspec:BadSpecsSize','The lower and upper limits of the shading area are badly specified');
end

if ~strcmp(region,'inside') && ~strcmp(region,'outside')
    error('FSDA:distribspec:BadRegion','The region can be either "inside" or "outside"');
end


% Optional parameters

options = struct('userColor', 'b');

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)

    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end

    % Check if all the specified optional arguments were present
    % in structure options
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
userColor = options.userColor;

%swap the specs if they are reversed
if specs(1) > specs(2)
    specs = fliplr(specs);
end
lb = specs(1);
ub = specs(2);
lbinf = isinf(lb);
ubinf = isinf(ub);

if lbinf && ubinf
    error('FSDA:distribspec:BadSpecsInfinite','Both limits are infinite');
end


%% plot the distribution

% initialise figure
nspecfig = figure;
nspecaxes = axes;
set(nspecaxes, 'Parent', nspecfig);
set(nspecaxes,'Nextplot','add');
hh    = plot(x,y,'k-','LineWidth',1.5);
xlims = get(nspecaxes,'Xlim');

% compute the endpoints of the spec limit lines and plot limit lines
pll =  [xlims(1);xlims(1)];
ypll = [0;eps];
if lbinf
    ll =  pll;
    yll = ypll;
else
    ll =  [lb; lb];
    yll = [0; pdf(pd, lb)];
end

pul =  [xlims(2);xlims(2)];
ypul = [eps;0];
if ubinf
    ul =  pul;
    yul = ypul;
else
    ul  = [ub; ub];
    yul = [pdf(pd, ub); 0];
end

% plot the shade area
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
    title({'User-defined PDF' ; char(userpdf) ; strpar ; strprob ; ['$lb =$ ' num2str(lb) ' -- ' '$ub =$ ' num2str(ub)]}, 'interpreter' , 'latex' , 'FontSize', 14);

elseif ~isempty(usercdf)
    % CASE 2.2.b: this is for a user-defined distribution: cdf
    if parnum>0
        strpar = ['with argmin(rmse) parameter ' , '$\hat{\theta}$' '=' num2str(parhat) ' '];
    else
        strpar = '';
    end
    title({'User-defined CDF' ; char(usercdf) ; strpar ; strprob ; ['$lb =$ ' num2str(lb) ' -- ' '$ub =$ ' num2str(ub)]}, 'interpreter' , 'latex' , 'FontSize', 14);

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
