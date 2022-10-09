function [p, h] = distribspec(pd, specs, region, varargin)
%distribspec plots a probability density function between specification limits
%
% <a href="matlab: docsearchFS('distribspec')">Link to the help function</a>
%
%    distribspec generalises the MATLAB function normspec.m for plotting a
%    selected probability density function by shading the portion inside
%    given limits.
%
% Required input arguments:
%
%    pd     : probability density. Object returned by makedist.
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
%    userColor :   The color of the shaded area.
%
% Output:
%
%    p:   The probability of the shaded area. Scalar in [0 1].
%
%    h:   A handle to the line objects.
%
%
% Optional Output:
%
%   xxxx     : The xxx.
%
% More About:
%
% zzzzzz
%
% See also: normspec, makedist
%
%
% References:
%
% Copyright 2008-2021.
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
    % A Gamma with parameter values a = 3 and b = 1, in [2 3].
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [2 3];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region);
%}

%{
    % A Gamma with parameter values a = 3 and b = 1, outside [2 3].
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
    [p, h] = distribspec(pd, specs, region, 'userColor','r');
%}

%{
    % A Gamma with parameter values a = 3 and b = 1, up to -1.
    pd = makedist('Gamma','a',3,'b',1);
    specs  = [-inf -1];
    region = 'inside';
    [p, h] = distribspec(pd, specs, region, 'userColor','r');
%}

%% Beginning of code

% Checking inputs and invalid arguments

if nargin<3
    region='inside';
end

if numel(specs) ~= 2 || ~isnumeric(specs)
    error('FSDA:distribspec:BadSpecsSize','The lower and upper limits of the shading area are badly specified');
end

if ~strcmp(region,'inside') && ~strcmp(region,'outside')
    error('FSDA:distribspec:BadRegion','The region can be either "inside" or "outside"');
end

if ~(isobject(pd) && isprop(pd,'DistributionName'))
    error('FSDA:distribspec:BadDistrib','Bad distribution object: use makedist.m to generate one');
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
    options.varargin{i}=varargin{i+1};
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
    error(message('FSDA:distribspec:BadSpecsInfinite'));
end


%% plot the distribution 

% set of values based on the given probability density 
prob = (0.0002:0.0004:0.9998)';
x = icdf(pd,prob);
y = pdf(pd, x);

% initialise figure
nspecfig = figure;
nspecaxes = axes;
set(nspecaxes, 'Parent', nspecfig);
set(nspecaxes,'Nextplot','add');
hh    = plot(x,y,'b-');
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
        xfill = [pll;  x(k1); ll          ; ul;          x(k2); pul  ];
        yfill = [ypll; y(k1); flipud(yll) ; flipud(yul); y(k2); ypul ];
        fill(xfill,yfill,userColor);
end

%% compute p
if strcmp(region,'outside')
    if lbinf
        p = cdf(pd,-ub); % P(t > ub)
        strprob = ['Probability greater than upper bound is ' num2str(p)];
    elseif ubinf
        p = cdf(pd,lb); % P(t < lb)
        strprob = ['Probability smaller than lower bound is ' num2str(p)];
    else
        p = cdf(pd,lb) + cdf(pd,-ub); % P(t < lb) + Pr(t > ub)
        strprob = ['Probability in the selected outside region is ' num2str(p)];
    end
else
    if lbinf
        p = cdf(pd,ub); % P(t < ub)
        strprob = ['Probability lower than upper bound is ' num2str(p)];
    elseif ubinf
        p = cdf(pd,-lb); % P(t > lb)
        strprob = ['Probability greater than lower bound is ' num2str(p)];
    else
        p = cdf(pd,ub) - cdf(pd,lb);  % P(lb < t < ub)
        strprob = ['Probability in the selected region is ' num2str(p)];
    end
end

%%  Add title and labels to the plot and return the handles

numpar = length(pd.ParameterNames);
strpar = [pd.DistributionName ' with '];
for np=1:numpar
    strpar = [strpar , char(pd.ParameterNames(np)) '=' num2str(pd.ParameterValues(np)) ' ']; %#ok<AGROW>
end

title({strpar ; strprob}, 'interpreter' , 'latex' , 'FontSize', 14);
xaxis = refline(0,0);
set(xaxis,'Color','k');
ylabel('Density value', 'interpreter' , 'latex' , 'FontSize', 12);
xlabel('Critical value', 'interpreter' , 'latex' , 'FontSize', 12);

if nargout > 1
    h = [hh; hh1];
end

end
