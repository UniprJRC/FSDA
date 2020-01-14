function [Wt,pretain,varargout] = wthin(X,varargin)
%wthin thins a uni/bi-dimensional dataset
%
%<a href="matlab: docsearchFS('wthin')">Link to the help page for this function</a>
%
%   Computes retention probabilities and bernoulli (0/1) weights on the
%   basis of data density estimate.
%
%  Required input arguments:
%
%   X :          Input data. Vector or 2-column matrix. The structure
%                contains the uni/bi-variate data to be thinned on the
%                basis of a probability density estimate.
%
%
%  Optional input arguments:
%
%   bandwidth :  bandwidth value. Scalar. The bandwidth used to estimate
%                the density. It can be estimated from the data using
%                function bwe.
%                Data Types - scalar
%                Example - bandwidth,0.35
%
%   retainby  :  retention method. String. The function used to retain the
%                observations. It can be:
%                - 'inverse' , i.e. (1 ./ pdfe) / max((1 ./ pdfe)))
%                - 'comp2one' (default),  i.e. 1 - pdfe/max(pdfe))
%                Data Types - char
%                Example - 'method','comp2one'
%
%
%  Output:
%
%   Wt :        vector of Bernoulli weights. Vector. Contains 1 for retained
%               units and 0 for thinned units.
%               Data Types - single | double.
%
%   pretain :   vector of retention probabilities. Vector. These are the
%               probabilities that each point in X will be retained,
%               estimated using a gaussian kernel using function ksdensity.
%               Data Types - single | double.
%
%  Optional Output:
%
%   Xt :        vector of retained units. Vector. It is X(Wt,:).
%               Data Types - single | double.
%
%  See also ksdensity, mvksdensity, bwe
%
%
% References:
%
% Bowman, A.W. and Azzalini, A. (1997), "Applied Smoothing
% Techniques for Data Analysis", Oxford University Press.
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('wthin')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%
% Examples:
%
%{
    % Univariate thinning.
    clear all; close all;
    % The dataset is bi-dimensional and contain two collinear groups with
    % regression structure. One group is dense, with 1000 units; the second
    % has 100 units. Thinning in done according to the density of the values
    % predicted by the OLS fit.
    x1 = randn(1000,1);
    x2 = 8 + randn(100,1);
    x = [x1 ; x2];
    y = 5*x + 0.9*randn(1100,1);
    b = [ones(1100,1) , x] \ y;
    yhat = [ones(1100,1) , x] * b;
    plot(x,y,'.',x,yhat);

    %x3 = 0.2 + 0.01*randn(1000,1);
    %y3 = 40 + 0.01*randn(1000,1);
    %plot(x,y,'.',x,yhat,'--',x3,y3,'.');

    % thinning over the predicted values
    %[Wt,pretain] = wthin([yhat ; y3], 'retainby','comp2one');
    [Wt,pretain] = wthin(yhat, 'retainby','comp2one');

    figure;
    plot(x(Wt,:),y(Wt,:),'k.',x(~Wt,:),y(~Wt,:),'r.');
    drawnow;
    axis manual;
    title('univariate thinning over predicted ols values')
    clickableMultiLegend(['Retained: ' num2str(sum(Wt))],['Thinned:   ' num2str(sum(~Wt))]);

%}

%{
    % Bi-dimensional thinning.
    % Same dataset, but thinning is done on the original bi-variate data.
   
    plot(x,y,'.');

    % thinning over the original bi-variate data
    [Wt2,pretain2] = wthin([x,y]);

    plot(x(Wt2,:),y(Wt2,:),'k.',x(~Wt2,:),y(~Wt2,:),'r.');
    drawnow;
    axis manual;
    title('bivariate thinning')
    clickableMultiLegend(['Retained: ' num2str(sum(Wt2))],['Thinned:   ' num2str(sum(~Wt2))]);
%}

%{
    % Use of 'retainby' option.
    % Since the thinning on the original bi-variate data with the default
    % retention method ('inverse') removes too many units, let's try with
    % the less conservative 'comp2one' option.
   
    plot(x,y,'.');

    % thinning over the original bi-variate data
    [Wt2,pretain2] = wthin([x,y], 'retainby','comp2one');

    plot(x(Wt2,:),y(Wt2,:),'k.',x(~Wt2,:),y(~Wt2,:),'r.');
    drawnow;
    axis manual
    clickableMultiLegend(['Retained: ' num2str(sum(Wt2))],['Thinned:   ' num2str(sum(~Wt2))]);
    title('"comp2one" thinning over the original bi-variate data');
    
%}

%{
    % Optional output Xt.
    % Same dataset, the retained data are also returned using varagout option.
   
    % thinning over the original bi-variate data
    [Wt2,pretain2,RetUnits] = wthin([x,y]);
    % disp(RetUnits)
%}

%{
   % thinning on the fishery dataset.
    load fishery;
    X=fishery.data;
    % some jittering is necessary because duplicated units are not treated
    % in tclustreg: this needs to be addressed
    X = X + 10^(-8) * abs(randn(677,2));

    % thinning over the original bi-variate data
    [Wt3,pretain3,RetUnits3] = wthin(X ,'retainby','comp2one');
    figure;
    plot(X(Wt3,1),X(Wt3,2),'k.',X(~Wt3,1),X(~Wt3,2),'rx');
    drawnow;
    axis manual
    clickableMultiLegend(['Retained: ' num2str(sum(Wt3))],['Thinned:   ' num2str(sum(~Wt3))]);
    title('"comp2one" thinning on the fishery dataset');

%}

%{
    % univariate thinning with less than 100 units.
    % As the first examp[le above, but with less than 100 units in the data.
    x1 = randn(85,1);
    x2 = 8 + randn(10,1);
    x = [x1 ; x2];
    y = 5*x + 0.9*randn(95,1);
    b = [ones(95,1) , x] \ y;
    yhat = [ones(95,1) , x] * b;
    plot(x,y,'.',x,yhat,'--');

    % thinning over the predicted values
    [Wt,pretain] = wthin(yhat, 'retainby','comp2one');

    plot(x(Wt,:),y(Wt,:),'k.',x(~Wt,:),y(~Wt,:),'r.');
    drawnow;
    axis manual
    title('univariate thinning over ols values predicted on a small dataset')
    clickableMultiLegend(['Retained: ' num2str(sum(Wt))],['Thinned:   ' num2str(sum(~Wt))]);

%}

%% options

% for reasons of performance options are checked only if necessary
if nargin > 1
    
    options     = struct('retainby','comp2one','bandwidth',0);
    UserOptions = varargin(1:2:length(varargin));
    if ~isempty(UserOptions) && (length(varargin) ~= 2*length(UserOptions))
        error('FSDA:kdebiv:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    if nargin>1
        % Write in structure 'options' the options chosen by the user
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
    end
    
    % retention method: it can be
    % 'comp2one' (i.e. 1 - pdfe/max(pdfe))
    % 'inverse'  (i.e. (1 ./ pdfe) / max((1 ./ pdfe)))
    retainby    = options.retainby;
    
    % the bandwidth used to estimate the density
    bandwidth   = options.bandwidth;
    
    if ~isempty(UserOptions)
        retainby_types  = {'inverse' , 'comp2one'};
        if  isempty(retainby) || ~(ischar(retainby) && max(strcmp(retainby,retainby_types)))
            retainby = 'inverse';
        end
        if  ~isvector(bandwidth)
            bandwidth = 0;
        end
        if bandwidth > 0
            % if the user has chosen a bandwidth, we may want to provide
            % a support too. For the moment we leave it unbounded, which is
            % the default of ksdensity.
            support = 'unbounded';%'positive ';
            %Remark: if we want to provide the support for the density
            %estimation, then the support should include the data values
            %interval. The quantity 'e' that exceeds the interval should
            %not be too small, otherwise (if bandwidth/2 is larger than
            %'e') the estimates of the pdf at the extreme data values risk
            %to be completely wrong. We have empirically observed this
            %effect using the following support values, with factor = 6 and
            %factor = 2:
            % factor = 2; %factor = 6;
            % minX = min(X); maxX = max(X); e = (maxX-minX)/10^(factor);
            % support = [ (min(X)-e) , (max(X)+e) ];
        end
    end
else
    bandwidth = 0;
    retainby  = 0;
end

%% Compute the density along the predicted values

[~,d] = size(X);

% The if statements below, which control the application of function
% ksdensity, are to optimize the time execution of the density estimation.
% We have empirically observed that:
% 1. ksdensity is faster if bandwidth is not provided (counter-intuitive,
%    but due to the effect of the options checks).
% 2. ksdensity is slower if it is evaluated at specified values provided in
%    the second input argument, say 'E' in ksdensity(X,E,'Support',support).
% 3. By default, when the evaluation vector E is not provided, ksdensity
%    evaluates the pdf on 100 points (if d=1), independently of the size of
%    the input (estimation) vector X. If the evaluation points are much
%    more than 100, then ksdensity becomes extremely slower. For this
%    reason, we estimate the pdf on the points in X, we interpolate the
%    pdf, and finally evaluate the interpolated function on the points in
%    E. In our specific case, where E is equal to X, we leave ksdensity to
%    compute the pdf on a default sample of 100 points taken or estimated
%    from X; then, we interpolate and evaluate the pdf on the full X
%    sample. For a sample of 1000 units, we reduce the time execution of
%    one order of magnitude (from ~14 to ~1.5 seconds).
if d > 1
    if ~verLessThan('matlab','9.0')
        % for the moment no optimization is done (points 2 and 3 not addressed).
        if bandwidth == 0
            % Remark: by default ksdensity estimates the bandwidt with Scott's rule.
            [pdfe,xout,u]  = ksdensity(X,X);
        else
            [pdfe,xout,u]  = ksdensity(X,X,'Support',support,'bandwidth',bandwidth);
            %[pdfe,xout,u]  = ksdensity(X,X,'Support','positive','BoundaryCorrection','reflection','bandwidth',bandwidth);
            %[pdfe,xout,u]  = ksdensity(X,X,'bandwidth',bandwidth);
        end
    else
        [pdfedef,xout1,u]  = kdebiv(X,'pdfmethod','fsda');
        if verLessThan('matlab', '8.1')
            Fpdfe = TriScatteredInterp(xout1(:,1),xout1(:,2),pdfedef); %#ok<DTRIINT>
        else
            Fpdfe = scatteredInterpolant(xout1(:,1),xout1(:,2),pdfedef);
        end
        pdfe  = Fpdfe(X(:,1),X(:,2));
        xout = X;
    end
else
    % This is the univariate case. We address points 2 and 3.
    % We do first an estimate of the pdf on the default 100 points adopted
    % by ksdensity.
    if bandwidth == 0
        [pdfedef,xout1,u]  = ksdensity(X);
    else
        [pdfedef,xout1,u]  = ksdensity(X,'Support',support,'bandwidth',bandwidth);
    end
    % Then we interpolate the 100 estimated pdf values on the data input
    % points. This is independent from the number of evaluation points,
    % because we need an estimate of the pdf exactely on the data input
    % points, not on other points chosen, e.g., as a linearly spaced vector
    % xq = linspace((min(X)-e) , (max(X)+e) , n);
    pdfe = interp1(xout1,pdfedef,X);
    xout = X;
end

varargout{2} = u;
varargout{3} = xout;

% replace the zero or negative (in case of rounding problems) sampling
% probability with a very small value
pdfe(pdfe<10^(-15))=10^(-15);
% convert the density values into the vector of retention probabilities;
% sampling probability should be inversely proportional to the density, but
% different functions are possible.
if retainby == 0
    % in this case the user has not provided optional arguments and accept
    % all defaults. For performance reasons, the 'switch' statement is
    % skipped and the default 'comp2one' function is applied.
    pretain = 1 - pdfe/max(pdfe);
    % pretain = (1 ./ pdfe) / max((1 ./ pdfe));
else
    switch retainby
        case 'comp2one'
            % complement to 1
            pretain = 1 - pdfe/max(pdfe);
        case 'inverse'
            % inverse
            pretain = (1 ./ pdfe) / max((1 ./ pdfe));
    end
end

% Thinning: Xt is the retained vector; Wt are the indices of the retained
% points in the original data X.
[Xt , Wt] =  rthin(X , pretain);
varargout{1} = Xt;

end

%FScategory:UTISTAT
