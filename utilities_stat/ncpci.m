function ci=ncpci(x,fType,df,varargin)
% non centrality parameter confidence interval (taken from effect_of_size_toolbox)
%
%<a href="matlab: docsearchFS('ncpci')">Link to the help function</a>
%
% This function creates (using iteration) a two-sided confidence intervals for the
% noncentrality parameter of a noncentral $\chi^2$, $F$ or
% $t$ distribution with degrees of freedom $df$, given an abscissa value $x$.
% This function has been taken from the MATLAB toolbox 'Measures of effect
% Size' by Harald Hentschke and  Maik C. Stüttgen.
% https://www.mathworks.com/matlabcentral/fileexchange/32398-hhentschke-measures-of-effect-size-toolbox
% and has been slighlty modified to be included into the FSDA
% toolbox.
%
%
%  Required input arguments:
%
%       x   :     Non centrality parameter for which confidence interval is
%                 needed.
%                 Scalar double.
%                 This is typically a value from $\chi^2$, $F$ or $t$
%                 distribution. Confidence interval based on  the
%                 noncentrality parameter of a noncentral
%                 distribution describes the degree of deviation from the
%                 null hypothesis. Its value is zero if the null hypothesis
%                 is true, and different from zero otherwise. For further
%                 details see the "More about" section of this document.
%     fType  :    distribution to use. Character.
%                 Character equal to 'X2', 'F' or 't' which specifies the
%                 reference distribution to use. X2 means $\chi^2$
%                 distribution, $F$ means Fisher F distribution and t means
%                 Student $t$ distribution.
%       df   :    degrees of freedom. Scalar or vector of length 2.
%                 Degrees of freedom of the reference distribution. It
%                 fType is 'X2' or 't' df must be a scalar. If, on the other
%                 hand, fType is 'F', then df must be a vector of length 2
%                 (which contains respectively the degrees of freedom of
%                 the numerator and of the denominator of the F
%                 distribution)
%
%
%  Optional input arguments:
%
% confLevel:     confidence levels to be used to
%               compute confidence intervals. Scalar.
%               The default value of conflev is 0.95  that
%               is 95 per cent confidence interval
%               is computed.
%               Example - 'confLevel',0.99
%               Data Types - double
%     prec   : tolerance for the iterative loop. Scalar.
%             Iteration will run until the estimated percentile is <=prec
%             away from the requested percentile.
%             The default value is 1e-6;
%                 Example - 'prec',1e-05
%                 Data Types - single | double
% doAnimate :   show graphically the iteration process.
%               Logical.
%               If doAnimate is true the the iteration process will be
%               graphically displayed in a figure window. The default value
%               of doAnimate is false.
%                 Example - 'doAnimate',false
%                 Data Types - logical
%
%  Output:
%
%        ci  :  confidence interval for the non centrality parameter.
%                 1-by-2 vector.
%                 Vector which contains the lower and upper confidence
%                 interval of the non centrality parameter.
%
%
%
% More About:
%
% This function is used in the FSDA toolbox in function corrNominal.m to
% find the confidence interval of Cramer's $V$ index. This index is a
% function of the non centrality parameter associated with the $\chi^2$
% index. Confidence intervals based on non central distributions depend on
% the "inversion confidence interval principle" (Stiegler and Fouladi 1997,
% pp. 237-239). The main idea is to use the observed value of a test
% statistic (that is input x in the language of this routine) to initiate a
% search for the lower and upper limits of to a $1-\alpha$ confidence
% interval for the non centrality parameter (the lower and upper bound of
% the confidence interval is given in output argument ci).
% The confidence interval of the non centrality parameter can then be
% converted into a confidence interval of an index which takes into account
% the sample size (which for example can be the Cramer's $V$ index) as long
% as the effect-size index (parameter) is a monotonic function of the non
% centrality parameter.
% See Smithson (2003) and "Measures of Effect Size" Toolbox for further
% details.
%
% See also: corrNominal, corrOrdinal, ncx2cdf, ncfcdf, nctcdf
%
% References:
%
% Hentschke, H. and Stüttgen, M. (2011), Comuputation of measures of effect
% size for neuroscience data sets, "European Journal of Neuroscience", Vol.
% 34, pp. 1887-1894.
% Smithson, M.J. (2003), "Confidence Intervals", Quantitative Applications in
% the Social Sciences Series, No. 140. Thousand Oaks, CA: Sage. [pp. 39-41]
% Hentschke, H. and Stüttgen, M. (2015), Measures of Effect Size Toolbox
% Version 1.4. [Code by Harald Hentschke (University of Tübingen) and
% Maik Stüttgen (University of Bochum)].
% Steiger, J.H., and Fouladi, R.T. (1997), Noncentrality interval
% estimation and the evaluation of statistical models. In Harlow, L.L.,
% Stanley, S., Mulaik, A. and Steiger, J.H., Eds., "What if there were no
% significance tests?", pp. 221-257. Mahwah, NJ: Lawrence Erlbaum.
%
%
% Acknowledgements:
%
% This function has been taken and adapted from the MATLAB toolbox 'Measures of effect
% Size' by Harald Hentschke and  Maik C. Stüttgen.
% https://www.mathworks.com/matlabcentral/fileexchange/32398-hhentschke-measures-of-effect-size-toolbox
%
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('ncpci')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%

% Examples:

%{
    %% ncpci with all the default values.
    % Suppose that in a contingency table of size 3-times-4 the value of%
    % the chi square test is 52. Suppose we want to compute a confidence
    % interval for the non centrality parameter of the chi^2 with
    % (3-1)(4-1)=8 degrees of freedom.
    ci=ncpci(52,'X2',8);
    disp('Confidence interval for the non centrality parameter')
    disp(ci)
%}

%{
    %% ncpci with option confint.
    % A 99 per cent confidence interval is requested.
    confint=0.99;
    ci=ncpci(52,'X2',8,'confLevel',confint);
    disp([ num2str(100*confint) ' per cent confidence interval for the non centrality parameter'])
    disp(ci)
%}

%{
    %% ncpci with option prec.
    % Increase the precision.
    prec=1e-12;
    ci=ncpci(52,'X2',8,'prec',prec);
    disp(['95 per cent confidence interval for the non centrality parameter'])
    disp(ci)
%}


%{
    % ncpci with option doAnimate.
    % set doAnimate to true.
    doAnimate=true;
    ci=ncpci(52,'X2',8,'doAnimate',true);
%}

%{
    % Confidence interval based on the F distibution.
    % See the animation which leads to convergence.
    ci=ncpci(52,'F',[8 3],'doAnimate',true);
%}


%% Beginning of code

% Defaults
prec=1e-6;
confLevel=.95;
doAnimate=false;

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    options=struct('prec',prec,...
        'confLevel',confLevel,'doAnimate',false);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:CorrNominal:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    % Write in structure 'options' the options chosen by the user
    if nargin > 2
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
    end
    prec=options.prec;
    confLevel=options.confLevel;
    doAnimate=options.doAnimate;
end

% Check that x is a positive number if fType is 'X2' or 'F'
isPosPdf=ismember(fType,{'X2','F'});
if isPosPdf && x<0
    error('FSDA:ncpci:WrongInputArg','input arg ''x'' is negative but must be positive for X2 and F distributions')
end

% Check that df is a vector of length 2 if fType is 'F'
isPosPdf=strcmp(fType,'F');
if isPosPdf && length(df)~=2
    error('FSDA:ncpci:WrongInputArg','F distribution has been specified and therefore input df must be a vector of lenght 2')
end

% convert df to cell for automatic expansion of parameters
df=num2cell(df);
% convert confidence level to alpha
alpha=1-confLevel;
% target p values
pTarget=[1-alpha/2  alpha/2];

% start index for outermost loop below, determining whether lower CI shall
% be computed or not
loopStartIx=1;
switch fType
    case 'X2'
        curPdf=@ncx2pdf;
        curCdf=@ncx2cdf;
        curInv=@chi2inv;
        % abscissa limits for plots (if doAnimate==true): first row for lower
        % CI, second row for upper CI
        abscissLim=[0 2*x;0 5*x];
        % check: if cdf of x with noncentrality parameter 0 is less than
        % 1-alpha/2 don't even start on the lower CI because the iteration will
        % not converge (that is, there is no lower CI for given values of x and
        % df)
        if chi2cdf(x,df{:})<1-alpha/2
            % lower CI cannot be constructed as it is too close to zero - set to
            % NaN
            ci=nan;
            loopStartIx=2;
        end
        
    case 'F'
        curPdf=@ncfpdf;
        curCdf=@ncfcdf;
        curInv=@finv;
        abscissLim=[0 2*x;0 5*x];
        % similar check as above
        if fcdf(x,df{:})<1-alpha/2
            % lower CI cannot be constructed as it is too close to zero - set to
            % NaN
            ci=nan;
            loopStartIx=2;
        end
        
    case 't'
        curPdf=@nctpdf;
        curCdf=@nctcdf;
        curInv=@tinv;
        abscissLim=x+[-4 2;-2 4]*sqrt(abs(x));
        
    otherwise
        error('FSDA:ncpci:InvalidArg','illegal distribution function specified');
end
%
% if prec>.001
%   warning('results will be inaccurate - set input parameter ''prec'' to a lower value');
% end

if doAnimate == true
    fh=figure;
    ph0=plot(x,0,'k^');
    hold on
    set(ph0,'markerfacecolor','k','markersize',6);
    ph=[];
    ti={'lower CI','upper CI'};
end

% loop twice: first lower ci (but see above), then upper ci
for iIx=loopStartIx:2
    % determine initial values: there are probably better ways of estimating
    % the limits of ncp for X2 and F pdfs than the guesses below (which work
    % best if the X2/F/t value is small)
    switch fType
        case 'X2'
            if iIx==1
                % lower CI
                ncp=x+curInv(pTarget(iIx),df{:});
            else
                % upper CI
                ncp=5*x;
            end
            
        case 'F'
            if iIx==1
                ncp=x+curInv(pTarget(iIx),df{:});
            else
                ncp=10*x;
            end
            
        case 't'
            % as a rough first approximation, assume that lower/upper limit of
            % ncp is close to corresponding percentiles of central pdfs
            if iIx==1
                ncp=x+curInv(pTarget(iIx),df{:});
            else
                ncp=x-curInv(pTarget(iIx),df{:});
            end
    end
    
    % interval of first estimates: guessed ncp enlarged by x/2 on either side
    ncp=ncp+abs(x)*[-.5 .5];
    % p values of current estimates
    p=curCdf(x,df{:},ncp);
    % deviations of p of current noncentral x pdfs from target p value
    deltaP=p-pTarget(iIx);
    nIter=1;
    if doAnimate == true
        ph=plotPdf(ncp,ph,curPdf,df,iIx,nIter,abscissLim,ti);
    end
    % while desired precision is not reached...
    while ~any(abs(deltaP)<=prec) && ~any(isnan(deltaP))     
        if all(deltaP>0)
            % shift interval to the right by one interval length
            ncp=[ncp(2) ncp(2)+abs(diff(ncp))];
        elseif all(deltaP<0)
            % shift left by one interval length
            ncp=[ncp(1)-abs(diff(ncp)) ncp(1)];
        else
            % halve interval around mean
            ncp=mean(ncp)+.25*abs(diff(ncp))*[-1 1];
        end
        % X2 and F distributions need an extra check: the lower ncp must be >=0
        if isPosPdf
            if ncp(1)<0
                ncp(1)=0;
            end
            % if both values of ncp are zero here the upper CI is zero, too, so
            % stop here
            if ~any(ncp)
                break
            end
        end
        % p values of current estimates
        p=curCdf(x,df{:},ncp);
        % deviations of p of current nc x pdfs from target
        deltaP=p-pTarget(iIx);
        nIter=nIter+1;
        if doAnimate == true
            ph=plotPdf(ncp,ph,curPdf,df,iIx,nIter,abscissLim,ti);
        end
    end
    % pick border which is closer to the target value
    [~,ix]=min(abs(deltaP));
    ci(iIx)=ncp(ix);
end

% close figure
if doAnimate == true
    pause(1)
    close(fh)
end
end

% ======================== inner function =================================
function ph=plotPdf(ncp,ph,pdfH,df,iIx,nIter,abscissLim,ti)
% ** function ph=plotPdf(ncp,ph,pdfH,df,iIx,nIter,abscissLim,ti)
% If doAnimate==true, plotPdf plots x (first input arg to ncpci) and
% noncentral pdfs with the noncentrality parameter estimates of each
% iteration step
abscissVal=linspace(abscissLim(iIx,1),abscissLim(iIx,2),200);
if isempty(ph)
    ph(1)=plot(abscissVal,pdfH(abscissVal,df{:},ncp(1)),'-');
    ph(2)=plot(abscissVal,pdfH(abscissVal,df{:},ncp(2)),'-');
    set(ph(1),'color',[.9 .3 .3]);
    set(ph(2),'color',[.3 .3 .9]);
else
    set(ph(1),'xdata',abscissVal,'ydata',pdfH(abscissVal,df{:},ncp(1)));
    set(ph(2),'xdata',abscissVal,'ydata',pdfH(abscissVal,df{:},ncp(2)));
end
title([ti{iIx} ', iteration # ' int2str(nIter)])
% supposedly, in animation mode we would like to be able to follow the
% iterative process with our eyes, so slow things down
drawnow
pause(0.1)
end
%FScategory:UTISTAT