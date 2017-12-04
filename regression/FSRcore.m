function [out]=FSRcore(INP,model,options)
%FSRcore scans the trajectory of mdr to check for exceedances
%
%
% Required input arguments:
%
%    INP    :   Structure containing monitoring of mdr and other quantities. Structure.
%               Structure containing the following fields.
%   INP.y   =   Response variable. Vector. Response variable, specified as
%               a vector of length n, where n is the number of
%               observations. Each entry in y is the response for the
%               corresponding row of X.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%  INP.X =      Predictor variables. Matrix. Matrix of explanatory
%               variables (also called 'regressors') of dimension n x p
%               where p denotes the number of explanatory variables
%               including the intercept.
%  INP.n =      Number of observations. Scalar. Number of rows of matrix X.
%  INP.p =      Number of predictor variables. Scalar. Number of columns of
%               matrix X.
%  INP.mdr =    Minimum deletion residual. Matrix.  n -init x 2 matrix
%               which contains the monitoring of minimum deletion residual
%               at each step of the forward search.
%               1st col = fwd search index (from init to n-1).
%               2nd col = minimum deletion residual.
%               Depending on the string 'model', mdr refers to OLS
%               mdr, GLS mdr or Bayes regression mdr.
%  INP.init =   Search initialization. Scalar.
%               It specifies the point where the user has started
%               monitoring mdr.
%  INP.Un  =    Units included in each step. Matrix.
%               (n-init) x 11 matrix which contains the unit(s) included
%               in the subset at each step of the fwd search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one. Un(1,2) for example contains
%               the unit included in step init+1. Un(end,2) contains the
%               units included in the final step of the search.
%               Un has 11 columns because we store up to 10 units
%               simultaneously in each step.
%   INP.bb=     Units included in each step. Matrix.
%               n-by-(n-init+1) or matrix n-by-r matrix which the units
%               belonging to the subset at each step of the forward search
%               or at selected steps.
%  INP.Bcoeff=  Estimated regression coefficients. Matrix.
%               (n-init+1) x (p+1) matrix containing the monitoring of
%               estimated beta coefficients in each step of the forward
%               search. The first column contains the fwd search index.
%               Depending on the string 'model', Bcoeff refers to OLS
%               coefficents, GLS coefficients or Bayes regression
%               coefficients.
% INP.Hetero=  Estimated coefficients in the skedastic equation. Matrix.
%               (n-init+1) x (r+1) matrix containing the monitoring of
%               estimated skedastic coefficients in each step of the forward
%               search. The first column contains the fwd search index.
%               This input is used just if strcmp(model,'H')
%  INP.S2 =   Estimate of $sigma^2$.  (n-init+1)-by-2 matrix containing the
%               monitoring of S2.  Depending on the string 'model', S2
%               refers to OLS, GLS or in the Baysian case it is the
%               posterior estimate of $\sigma^2$.
%   INP.Z =     Predictor variables in the regression equation (necessary
%               input just if model='H'). Matrix.
%               n x r matrix or vector of length r.
%               If Z is a n x r matrix it contains the r variables which
%               form the scedastic function as follows (if input option art==1)
%               \[
%               \omega_i = 1 + exp(\gamma_0 + \gamma_1 Z(i,1) + ...+ \gamma_{r} Z(i,r))
%               \]
%               If Z is a vector of length r it contains the indexes of the
%               columns of matrix X which form the scedastic function as
%               follows
%               \[
%               \omega_i = 1 +  exp(\gamma_0 + \gamma_1 X(i,Z(1)) + ...+
%               \gamma_{r} X(i,Z(r)))
%               \]
%
%   INP.beta0 = Prior mean of $\beta$ (necessary
%               input just if model='B'). p-times-1 vector.
%   INP.R     = Matrix associated with covariance matrix of $\beta$ (necessary
%               input just if model='B'). p-times-p
%               positive definite matrix.
%               It can be interpreted as X0'X0 where X0 is a n0 x p
%               matrix coming from previous experiments (assuming that the
%               intercept is included in the model)
%
%               The prior distribution of $\tau_0$ is a gamma distribution with
%               parameters $a_0$ and $b_0$, that is
%
%                \[     p(\tau_0) \propto \tau^{a_0-1} \exp (-b_0 \tau)
%                       \qquad   E(\tau_0)= a_0/b_0               \]
%
%   INP.tau0 =  Prior estimate of tau (necessary
%               input just if model='B'). Scalar. Prior estimate of $\tau=1/ \sigma^2 =a_0/b_0$.
%     INP.n0 =  Number of previous experiments (necessary
%               input just if model='B'). Scalar. Sometimes it helps
%               to think of the prior information as coming from n0
%               previous experiments. Therefore we assume that matrix X0
%               (which defines R), was made up of n0 observations.
%                 Data Types - struct
%  model :      type of regression model. Character.
%               Possible values are '' (default) | 'H' | 'B'.
%               '' stands for linear regression;
%               'H' stands for heteroskedastic regression;
%               'B' stands for Bayesian regression.
%               This input is used to reconstruct the units belonging to
%               subset at step n-decl where decl is the number of units
%               declared as outliers. More precisely, if n>5000 matrix BB
%               just contains the units belonging to subset in selected
%               steps, therefore in order to find the units inside subset at
%               step n-decl, FSRcore calls:
%               routine FSRbsb.m in presence of linear regression;
%               routine FSRHbsb.m in presence of heteroskedastic regression;
%               routine FSRBbsb.m in presence of Bayesian regression;
%                 Data Types - char
%    options:   Additional options. Stucture. Structure containing optional
%               parameters which are passed to directly through functions
%               FSR.m, FSRH.m or FSRB.m.
%                 Data Types - struct
%
% Optional input arguments:
%
% Output:
%
%     out :     A structure containing the following fields
% out.ListOut  =  k x 1 vector containing the list of the units declared as
%                 outliers or NaN if the sample is homogeneous.
% out.outliers =  out.ListOut. This field is added for homogeneity with the
%                 other robust estimators.
% out.beta   =  p-by-1 vector containing the estimated regression
%               parameter in step n-k. Depending on the string 'model',
%               beta refers to OLS coefficents, GLS coefficients or Bayes
%               regression coefficients.
% out.scale   = estimate of the scale. Depending on the string 'model',
%               beta refers to OLS coefficents, GLS coefficients or it is
%               the inverse of the posterior estimate of the square root of tau.
% out.mdr    =  (n-init) x 2 matrix
%               1st col = fwd search index
%               2nd col = value of minimum deletion residual in each step
%               of the fwd search. Depending on the string 'model',
%               mdr is found using linear regression, heteroskedastic
%               regression or Bayes regression.
% out.Un     =  (n-init) x 11 matrix which contains the unit(s) included
%               in the subset at each step of the fwd search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one.
%               Un(1,2) for example contains the unit included in step
%               init+1.
%               Un(end,2) contains the units included in the final step
%               of the search.
% out.nout    = 2 x 5 matrix containing the number of times mdr went out
%               of particular quantiles.
%               First row contains quantiles 1 99 99.9 99.99 99.999.
%               Second row contains the frequency distribution.
%
% More About:
%
% The rules for declaring units as outliers are the same for standard
% regression, heteroskedastic regression and Bayesian regression. Therefore
% this function is called by:
% FSR.m  = outlier detection procedure for linear regression;
% FSRB.m = outlier detection procedure in Bayesian linear regression;
% FSRH.m = outlier detection procedure for heteroskedastic models;
% If ndecl units are declared as outliers, it is necessary to find the
% units forming subset at step n-decl. If n<=5000 input matrix INP.bb
% contains the storing of the units belonging to subset in all steps, else
% if INP.bb does not contain the units in step n-decl procedure calls
% routine FSRbsb.m or FSRHbsb.m or FSRBbsb.m.
%
% See also: FSR.m, FSRH.m, FSRB.m
%
% References:
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
% Last modified 31-05-2016

%% Beginning of code

y=INP.y;
X=INP.X;
n=INP.n;
p=INP.p;
mdr=INP.mdr;
init=INP.init;
Un=INP.Un;
bb=INP.bb;
Bcoeff=INP.Bcoeff;

S2=INP.S2;

% intcolumn = the index of the first constant column found in X, or empty.
% Used here to check if X includes the constant term for the intercept.
% The variable 'intercept' will be used later for plotting.
intcolumn = find(max(X,[],1)-min(X,[],1) == 0,1);
if any(intcolumn) && p>1
    intercept=1;
else
    intercept=0;
end


plo=options.plots;
labeladd=options.labeladd;
bivarfit=options.bivarfit;
multivarfit=options.multivarfit;
xlimx=options.xlim;
ylimy=options.ylim;
msg=options.msg;

bonflev=options.bonflev;
seq=1:n;

% correction in case of Bayesian model to account for number of
% (fictitious) observations in the prior
if strcmp(model,'B')
    nori = n;
    n    = nori + INP.n0;
end

if ~isempty(bonflev)
    if bonflev<1
        [gbonf] = FSRbonfbound(n,p,'prob',bonflev,'init',init);
        bonfthresh=gbonf;
    else
        bonfthresh=bonflev*ones(n-init,1);
    end
    
    % correction in case of Bayesian model to account for prior
    % observations
    if strcmp(model,'B')
        bonfthresh(1:INP.n0,:)=[];
        bonfthresh(:,1) = bonfthresh(:,1) - INP.n0;
        n = nori;
    end
    
else
    
    exact=1;
    % lowexceed=0 means than outlier detection is just based on upper
    % exceedances
    lowexceed=0;
    
    % Compute theoretical envelopes based on all observations
    quant=[0.99;0.999;0.9999;0.99999;0.01;0.5;0.00001];
    % Compute theoretical envelopes for minimum deletion residual based on all
    % the observations for the above quantiles.
    [gmin] = FSRenvmdr(n,p,'prob',quant,'init',init,'exact',exact);
    % gmin = the matrix which contains envelopes based on all observations.
    % 1st col of gmin = fwd search index
    % 2nd col of gmin = 99% envelope
    % 3rd col of gmin = 99.9% envelope
    % 4th col of gmin = 99.99% envelope
    % 5th col of gmin = 99.999% envelope
    % 6th col of gmin = 1% envelope
    % 7th col of gmin = 50% envelope
    
    % correction in case of Bayesian model to account for number of
    % (fictitious) observations in the prior
    if strcmp(model,'B')
        gmin(1:INP.n0,:)=[];
        gmin(:,1) = gmin(:,1) - INP.n0;
        n = nori;
    end
    
    % Thus, set the columns of gmin where the theoretical quantiles are located.
    [c99 , c999 , c9999 , c99999 , c001 , c50] = deal(2,3,4,5,6,7);
    
    
    bool=mdr(:,1)>=init;
    mdr=mdr(bool,:);
    gmin=gmin(gmin(:,1)>=mdr(1,1),:);
    
    
    % Store in nout the number of times the observed mdr (d_min) lies above:
    [out99 , out999 , out9999 , out99999 , out001] = deal( ...
        mdr(mdr(:,2)>gmin(:,c99),:) , ...       % the 99% envelope
        mdr(mdr(:,2)>gmin(:,c999),:) , ...      % the 99.9% envelope
        mdr(mdr(:,2)>gmin(:,c9999),:) , ...     % the 99.99% envelope
        mdr(mdr(:,2)>gmin(:,c99999),:) , ...    % the 99.999% envelope
        mdr(mdr(:,2)<gmin(:,c001),:) );         % the 1% envelope
    
    nout = [[1 99 999 9999 99999]; ...
        [size(out001,1) size(out99,1) size(out999,1) size(out9999,1) size(out99999,1)]];
    
    % NoFalseSig = boolean linked to the fact that the signal is good or not
    NoFalseSig=0;
    
    % NoFalseSig is set to 1 if the condition for an INCONTROVERTIBLE SIGNAL is
    % fulfilled.
    n9999 = nout(2,nout(1,:)==9999);
    if (n9999>=10)
        NoFalseSig=1;
        if msg
            disp('Observed curve of r_min is at least 10 times greater than 99.99% envelope'); % exact number is int2str(n9999)
            disp('--------------------------------------------------');
        end
    end
    
    % Divide central part from final part of the search
    istep = n-floor(13*sqrt(n/200));
end


%% Part 1. Signal detection and validation
nmdr=size(mdr,1);
if nmdr<4  && isempty(bonflev)
    error('FSDA:FSRcore:TooSmallRationp','Ratio n/p too small; modify init (i.e. decrease initial subset size) or use option bonflev (i.e. ''bonflev'', 0.99)')
end
signal=0;
sto=0;
extram3='';
extram2='';
strplot='';
resup=2;
if msg
    disp('-------------------------')
    disp('Signal detection loop');
end

%% Stage 1a: signal detection
% Signal dection is based on monitoring consecutive triplets or single
% extreme values

% Check if signal must be based on consecutive exceedances of envelope
% of mdr or on exceedance of global Bonferroni level
if isempty(bonflev)
    % Signal detection loop
    for i=3:nmdr
        
        if i<istep-init+1 % CENTRAL PART OF THE SEARCH
            % Extreme triplet or an extreme single value
            % Three consecutive values of d_min above the 99.99% threshold or 1
            % above 99.999% envelope
            if ((mdr(i,2)>gmin(i,c9999) && mdr(i+1,2)>gmin(i+1,c9999) && mdr(i-1,2)>gmin(i-1,c9999)) || mdr(i,2)>gmin(end,c99) || mdr(i,2)>gmin(i,c99999))
                if msg
                    disp(['Tentative signal in central part of the search: step m=' int2str(mdr(i,1)) ' because']);
                end
                if (mdr(i,2)>gmin(i,c9999) && mdr(i+1,2)>gmin(i+1,c9999) && mdr(i-1,2)>gmin(i-1,c9999))
                    if msg
                        disp(['rmin('  int2str(mdr(i,1)) ',' int2str(n) ')>99.99% and rmin(' int2str(mdr(i-1,1)) ',' int2str(n) ')>99.99% and rmin(' int2str(mdr(i+1,1)) ',' int2str(n) ')>99.99%']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99.99\%$ and $r_{min}(' int2str(mdr(i-1,1)) ',' int2str(n) ')>99.99\%$ and $r_{min}(' int2str(mdr(i+1,1)) ',' int2str(n) ')>99.99\%$'];
                    mdrsel=mdr(i-1:i+1,1:2);
                end
                
                if (mdr(i,2)>gmin(i,c99999))
                    if msg
                        disp(['rmin(' int2str(mdr(i,1)) ',' int2str(n) ')>99.999%']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99.999\%$'];
                    mdrsel=mdr(i-1:i+1,1:2);
                end
                
                if (mdr(i,2)>gmin(end,c99))
                    if msg
                        disp(['rmin(' int2str(mdr(i,1)) ',' int2str(n) ')>99% at final step: Bonferroni signal in the central part of the search.']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99\%$ at final step (Bonferroni signal)'];
                    mdrsel=mdr(i,1:2);
                    NoFalseSig=1; % i.e., no need of further validation
                end
                
                '------------------------------------------------';
                
                signal=1;
                
            elseif (mdr(i,2)<gmin(i,end)) && lowexceed==1  % && mdr(i,1)>round(n/2); % exceedance of the lower band
                if msg
                    disp(['rmin(' int2str(mdr(i,1)) ',' int2str(n) ')<0.00001% Bonferroni signal in the central part of the search.']);
                end
                strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')<0.00001\%$  (Bonferroni signal)'];
                mdrsel=mdr(i,1:2);
                NoFalseSig=1; % i.e., no need of further validation
                signal=2;
                mdag=mdr(i,1);
            else
            end
        elseif i<size(mdr,1)-1 % FINAL PART OF THE SEARCH
            % Extreme couple adjacent to an exceedance
            % Two consecutive values of mdr above the 99.99% envelope and 1 above 99%
            %             if ( (mdr(i,2)   > gmin(i,c999)   && mdr(i+1,2) > gmin(i+1,c999) && mdr(i-1,2) > gmin(i-1,c99)) || ...
            %                  (mdr(i-1,2) > gmin(i-1,c999) && mdr(i,2)   > gmin(i,c999)   && mdr(i+1,2) > gmin(i+1,c99)) || ...
            %                   mdr(i,2)   > gmin(end,c99)  || ...
            %                   mdr(i,2)   > gmin(i,c99999) )
            
            condition1 = (mdr(i,2)   > gmin(i,c999)   && mdr(i+1,2) > gmin(i+1,c999) && mdr(i-1,2) > gmin(i-1,c99)) ;
            condition2 = (mdr(i-1,2) > gmin(i-1,c999) && mdr(i,2)   > gmin(i,c999)   && mdr(i+1,2) > gmin(i+1,c99)) ;
            condition3 = mdr(i,2)   > gmin(end,c99)  ;
            condition4 = mdr(i,2)   > gmin(i,c99999) ;
            if ( condition1 || condition2 || condition3  || condition4 )
                
                if msg
                    disp(['Signal in final part of the search: step ' num2str(mdr(i,1)) ' because']);
                end
                
                if condition1
                    if msg
                        disp(['rmin('  int2str(mdr(i,1)) ',' int2str(n) ')>99.9% and rmin('  int2str(mdr(i+1,1)) ',' int2str(n) ')>99.9% and rmin('  int2str(mdr(i-1,1)) ',' int2str(n) ')>99%']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99.9\%$ and $r_{min}(' int2str(mdr(i-1,1)) ',' int2str(n) ')>99\%$ and $r_{min}(' int2str(mdr(i+1,1)) ',' int2str(n) ')>99.9\%$'];
                    mdrsel=mdr(i-1:i+1,1:2);
                end
                
                if condition2
                    if msg
                        disp(['rmin('  int2str(mdr(i-1,1)) ',' int2str(n) ')>99.9% and rmin('  int2str(mdr(i,1)) ',' int2str(n) ')>99.9% and rmin('  int2str(mdr(i+1,1)) ',' int2str(n) ')>99%']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99.9\%$ and $r_{min}(' int2str(mdr(i-1,1)) ',' int2str(n) ')>99.9\%$ and $r_{min}(' int2str(mdr(i+1,1)) ',' int2str(n) ')>99\%$'];
                    mdrsel=mdr(i-1:i+1,1:2);
                end
                
                if condition3
                    if msg
                        disp(['rmin('  int2str(mdr(i,1)) ',' int2str(n) ')>99% at final step: Bonferroni signal in the final part of the search.']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99\%$ at final step (Bonferroni signal)'];
                    mdrsel=mdr(i:i,1:2);
                end
                
                % Extreme single value above the upper threshold
                if condition4
                    if msg
                        disp(['rmin('  int2str(mdr(i,1)) ',' int2str(n) ')>99.999%']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99.999\%$'];
                    mdrsel=mdr(i:i,1:2);
                end
                
                '------------------------------------------------';
                % Signal is always considered true if it takes place in the
                % final part of the search
                NoFalseSig=1;
                signal=1;
            end
        elseif (mdr(i,2)>gmin(i,c999) || mdr(i,2)>gmin(end,c99)) && i==size(mdr,1)-1
            % potential couple of outliers
            signal=1;
            if msg
                disp('Signal is in penultimate step of the search');
            end
            
            if (mdr(i,2)>gmin(i,c999))
                if msg
                    disp(['rmin(' int2str(mdr(i,1)) ',' int2str(n) ')>99.9%']);
                end
                strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99.9\%$'];
            end
            
            if (mdr(i,2)>gmin(end,c99))
                if msg
                    disp(['rmin('  int2str(mdr(i,1)) ',' int2str(n) ')>99% at final step: Bonferroni signal in the final part of the search.']);
                end
                strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99\%$ at final step (Bonferroni signal)'];
            end
            mdrsel=mdr(i:i,1:2);
        elseif  mdr(i,2)>gmin(i,c99) && i==size(mdr,1)
            % a single outlier
            signal=1;
            if msg
                disp('Signal is in final step of the search');
            end
            strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99\%$ at final step'];
            mdrsel=mdr(i:i,1:2);
            
        elseif mdr(i,2)<gmin(i,end) && lowexceed==1 % exceedance of the lower extreme envelope
            signal=2;
            if msg
                disp(['rmin('  int2str(mdr(i,1)) ',' int2str(n) ')<0.0001%']);
            end
            strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')<0.0001\%$'];
            mdrsel=mdr(i:i,1:2);
            mdag=mdr(i,1);
        end
        
        %% Stage 1b: signal validation
        if signal==1
            if msg
                disp('-------------------')
                disp('Signal validation exceedance of upper envelopes');
            end
            % mdag is $m^\dagger$
            mdag=mdr(i,1);
            
            if mdr(i,1)<n-2
                % Check if the signal is incontrovertible
                % Incontrovertible signal = 3 consecutive values of d_min >
                % 99.999% threshold
                if mdr(i,2)>gmin(i,c99999) && mdr(i-1,2)>gmin(i-1,c99999) &&  mdr(i+1,2)>gmin(i+1,c99999)
                    if msg
                        disp(['3 consecutive values of r_min greater than 99.999% envelope in step mdag= ' int2str(mdr(i,1))]);
                    end
                    NoFalseSig=1;
                    extram3='Extreme signal';
                end
            else
                NoFalseSig=1;
            end
            
            % if the following statement is true, observed curve of r_min is
            % above 99.99% and later is below 1%: peak followed by dip
            if size(mdr,1)>mdag-mdr(1,1)+31
                if sum(mdr(i+1:i+31,2)<gmin(i+1:i+31,c001))>=2
                    NoFalseSig=1;  % Peak followed by dip
                    extram2='Peak followed by dip (d_min is above 99.99% threshold and in the sucessive steps goes below 1% envelope';
                end
            else
                if sum(mdr(i+1:end,2) < gmin(i+1:end,c001))>=2
                    NoFalseSig=1;  %Peak followed by dip in the final part of the search';
                    extram2='Peak followed by dip (r_min is above 99.99% threshold and in the sucessive steps goes below 1% envelope)';
                end
            end
            
            % if at this point NoFalseSig==0 it means that:
            % 1) n9999<10
            % 2) signal tool place in the central part of the search
            % 3) signal was not incontrovertible
            % 4) there was not a peak followed by dip
            if NoFalseSig==0
                % Compute the final value of the envelope based on
                % mdr(i+1,1)=mdagger+1 observations
                %[gval]=FSRenvmdr(mdag+1,p,'prob',0.01,'m0',mdag);
                [gval]=FSRenvmdr(mdag+1,p,'prob',0.01,'init',mdag);
                if mdr(i,2)<gval(1,2)
                    if msg
                        disp('false signal in step');
                        disp(['mdag='  int2str(mdag)]);
                    end
                    % increase mdag of the search by one unit
                    mdag=0;
                else
                    NoFalseSig=1;
                end
            end
            
            % If the signal has been validated get out of the signal detection
            % loop and move to stage 2: superimposition of the envelopes
            if (NoFalseSig==1)
                if msg
                    disp('Validated signal');
                end
                break
            end
        elseif signal==2
            if msg
                disp('-------------------')
                disp('Signal validation exceedance of extreme lower envelope');
            end
            break
        end
    end
else
    for i=1:nmdr
        % Outlier detection based on Bonferroni threshold
        if (mdr(i,2)>bonfthresh(i,end)) % TODO TODO && mdr(i,1)>floor(0.5*n);
            if msg
                disp(['mdr(' int2str(mdr(i,1)) ',' int2str(n) ')>99% Bonferroni level']);
            end
            strplot=['$mdr(' int2str(mdr(i,1)) ',' int2str(n) ')>99\%$ (Bonferroni level)'];
            mdrsel=mdr(i:i,1:2);
            
            signal=1;
            break
        end
    end
end

%% Create figure containing mdr + envelopes based on all the observations.
if plo==1 || plo ==2
    % get screen size [left, bottom, width, height]
    scrsz = get(0,'ScreenSize');
    
    figure1 = figure('Position',[1 scrsz(4)/2.5 scrsz(3)/3 scrsz(4)/2],'PaperSize',[20.98 29.68],'Name','Envelopes based on all the observations');
    axes1 = axes('Parent',figure1);
    
    % Create a pan-handle for the figure ...
    hpan_figure1 = pan(figure1);
    % ... and listen to pan events using callback functions
    set(hpan_figure1,'ActionPreCallback',@fprecallback);
    set(hpan_figure1,'ActionPostCallback',@fpostcallback);
    
    %     % Create a zoom-handle for the figure ...
    %     hzoom_figure1 = zoom(figure1);
    %     % ... and listen to pan events using callback functions
    %     set(hzoom_figure1,'ActionPreCallback',@figure1_precallback);
    %     set(hzoom_figure1,'ActionPostCallback',@figure1_postcallback_zoom);
    
    if isempty(xlimx)
        xl1=init-3; xl2=mdr(end,1);
    else
        xl1=xlimx(1);
        xl2=xlimx(2);
    end
    
    if isempty(ylimy)
        if isempty(bonflev)
            yl1=min([gmin(:,c001);mdr(:,2)]);
            yl2=max([gmin(:,c999);mdr(:,2)]);
        else
            yl1=min([bonfthresh(:,2);mdr(:,2)]);
            yl2=max([bonfthresh(:,2);mdr(:,2)]);
        end
    else
        yl1=ylimy(1);
        yl2=ylimy(2);
    end
    
    xlim([xl1 xl2]);
    ylim([yl1 yl2]);
    
    
    kx=0; ky=0;
    
    box('on'); hold('all');
    
    plot(mdr(:,1),mdr(:,2));
    
    
    if isempty(bonflev)
        
        % 1%, 99%, 99.9% envelopes based on all the observations
        line(gmin(:,1),gmin(:,[c001 c99 c999]),'Parent',axes1,'LineWidth',2,'LineStyle','--','Color',[0 0 1]);
        % 99.99% and 99.999% envelopes based on all the observations
        line(gmin(:,1),gmin(:,[c9999 c99999]),'Parent',axes1,'LineWidth',2,'LineStyle','--','Color',[1 0 0]);
        if signal==2
            line(gmin(:,1),gmin(:,end),'Parent',axes1,'LineWidth',2,'LineStyle','--','Color',[1 0 0]);
        end
        
        % 50% envelope based on all the observations
        line(gmin(:,1),gmin(:,c50),'Parent',axes1,'LineWidth',2,'LineStyle','--','Color',[1 0.69 0.39]);
        
        % Property-value pairs which are common to all quantile-labels
        PrVaCell{1,1} = 'HorizontalAlignment'; PrVaCell{2,1} = 'center';
        PrVaCell{1,2} = 'EdgeColor'; PrVaCell{2,2} = 'none';
        PrVaCell{1,3} = 'BackgroundColor'; PrVaCell{2,3} = 'none';
        PrVaCell{1,4} = 'FitBoxToText'; PrVaCell{2,4} = 'off';
        PrVaCell{1,5} = 'Tag'; PrVaCell{2,5} = 'quantile_label';
        
        % Create textbox with 1% label
        [figx, figy] = dsxy2figxy(gca, init, gmin(1,c001));
        if figy>=0 && figy<=1 && figx>=0 && figx<=1
            annotation(figure1,'textbox',[figx figy kx ky],...
                'String',{'1%'},...
                'UserData',[gmin(:,1) gmin(:,c001)],...
                PrVaCell{:});
        end
        
        % Create textbox with 99% label
        [figx, figy] = dsxy2figxy(gca, init, gmin(1,c99));
        
        if figy>=0 && figy<=1 && figx>=0 && figx<=1
            annotation(figure1,'textbox',[figx figy kx ky],...
                'String',{'99%'},...
                'UserData',[gmin(:,1) gmin(:,c99)],...
                PrVaCell{:});
        end
        
        % Create textbox with 50% label
        [figx, figy] = dsxy2figxy(gca, init, gmin(1,c50));
        
        if figy>=0 && figy<=1 && figx>=0 && figx<=1
            annotation(figure1,'textbox',[figx figy kx ky],...
                'String',{'50%'},...
                'UserData',[gmin(:,1) gmin(:,c50)],...
                PrVaCell{:});
        end
        
        if signal==1
            % Create textbox with 99.9% label
            [figx, figy] = dsxy2figxy(gca, init, gmin(1,c999));
            if figy>=0 && figy<=1 && figx>=0 && figx<=1
                annotation(figure1,'textbox',[figx figy kx ky],...
                    'String',{'99.9%'},...
                    'UserData',[gmin(:,1) gmin(:,c999)],...
                    PrVaCell{:});
            end
            
            % Create textbox with 99.99% label
            [figx, figy] = dsxy2figxy(gca, init, gmin(1,c9999));
            if figy<=1 && figy>=0 && figx>=0 && figx<=1
                annotation(figure1,'textbox',[figx figy kx ky],...
                    'String',{'99.99%'},...
                    'UserData',[gmin(:,1) gmin(:,c9999)],...
                    PrVaCell{:});
            end
            
            if gmin(1,c99999)<=yl2
                % Create textbox with 99.999% label
                [figx, figy] = dsxy2figxy(gca, init, gmin(1,c99999));
                if figy<=1 && figy>=0 && figx>=0 && figx<=1
                    annotation(figure1,'textbox',[figx figy kx ky],...
                        'String',{'99.999%'},...
                        'UserData',[gmin(:,1) gmin(:,c99999)],...
                        PrVaCell{:});
                end
            end
        else
        end
        
        % Add string which informs about the step where signal took place
        if signal==1 || signal==2
            strsig=['Signal is in step $m=' int2str(mdag) '$ because'];
            stem(mdrsel(:,1),mdrsel(:,2),'LineWidth',1,...
                'Color',[0.4784 0.06275 0.8941], 'DisplayName','Signal');
        else
            strsig='No signal during the search';
        end
        
        % Property-value pairs which are common to the next latex annotations
        PrVaCell{1,1} = 'Interpreter'; PrVaCell{2,1} = 'latex';
        PrVaCell{1,2} = 'HorizontalAlignment'; PrVaCell{2,2} = 'center';
        PrVaCell{1,3} = 'FitBoxToText'; PrVaCell{2,3} = 'on';
        PrVaCell{1,4} = 'EdgeColor'; PrVaCell{2,4} = 'none';
        PrVaCell{1,5} = 'BackgroundColor'; PrVaCell{2,5} = 'none';
        
        
        % latex annotations informing that the envelopes are based on
        % all the observations
        strmin=['$r_{min}(m,' int2str(n) ')$. '];
        annotation(figure1,'textbox',[0.5 0.9 kx ky],'String',{[strmin strsig]},...
            PrVaCell{:});
        
        annotation(figure1,'textbox',[0.5 0.85 kx ky],'String',{strplot},...
            PrVaCell{:});
        
        if istep<=xl2
            % Add vertical line which divides central part from final part of the
            % search
            [figx, figy] = dsxy2figxy(gca, istep, yl1);
            [figx, figy2] = dsxy2figxy(gca, istep, yl2);
            if figy2>=1
                figy2=1;
            else
            end
            if figx>=0
                annotation(figure1,'line',[figx figx],[figy figy2],...
                    'UserData',[istep yl1 yl2],...
                    'Tag','FinalPartLine');
            end
        end
    else
        % Superimpose Bonferroni line to the plot
        line(bonfthresh(:,1),bonfthresh(:,end),'Parent',axes1,'LineWidth',2,'LineStyle','--','Color',[0 0 1]);
        
        % Property-value pairs which are common to the next latex annotations
        PrVaCell=cell(2,5);
        PrVaCell{1,1} = 'Interpreter'; PrVaCell{2,1} = 'latex';
        PrVaCell{1,2} = 'HorizontalAlignment'; PrVaCell{2,2} = 'center';
        PrVaCell{1,3} = 'FitBoxToText'; PrVaCell{2,3} = 'on';
        PrVaCell{1,4} = 'EdgeColor'; PrVaCell{2,4} = 'none';
        PrVaCell{1,5} = 'BackgroundColor'; PrVaCell{2,5} = 'none';
        
        
        
        if size(bonfthresh,2)>1
            % latex annotations informing that the envelopes are based on
            % all the observations
            strmin='Exceedance based on Bonferroni threshold';
            annotation(figure1,'textbox',[0.5 0.9 kx ky],'String',strmin,...
                PrVaCell{:});
            msgth=['$r_{min}(' num2str(mdr(i,1)) ',' int2str(n) ')>' num2str(100*bonflev) '$\% envelope'];
            annotation(figure1,'textbox',[0.5 0.8 kx ky],'String',msgth,PrVaCell{:});
        else
            strmin='Exceedance based on user supplied threshold';
            annotation(figure1,'textbox',[0.5 0.9 kx ky],'String',strmin,...
                PrVaCell{:});
            msgth=['$r_{min}(' num2str(mdr(i,1)) ',' int2str(n) ')>$' num2str(bonflev)];
            annotation(figure1,'textbox',[0.5 0.8 kx ky],'String',msgth,PrVaCell{:});
        end
    end
    if signal==1
        stem(mdr(i,1),mdr(i,2),'LineWidth',1,...
            'Color',[0.4784 0.06275 0.8941], 'DisplayName','Signal');
    end
end

%% Part 2: envelope resuperimposition
% if a validated signal took place, superimposition of the envelopes starts
% from m^\dagger-1

if signal==1 || signal==2
    if isempty(bonflev) && signal==1
        
        if msg
            disp('-------------------------------');
            disp(['Start resuperimposing envelopes from step m=' int2str(mdag-1)]);
        end
        
        if plo==2
            % jwind is associated with subplot window number
            % nr is the number of row panes in the plot
            % nc is the number of columns panes in the plot
            jwind=1;
            nc=2;
            if mdr(i,1)>=n-2
                nr=1;
            else
                nr=2;
            end
            
            figure2 = figure('PaperSize',[20.98 29.68],'Name','Resuperimposed envelopes #1');
            % Create axes
            axes('Parent',figure2);
        end
        
        % First resuperimposed envelope is based on mdag-1 observations
        % Notice that mdr(i,1) = m dagger
        for tr=(mdag-1):(n)
            % Compute theoretical envelopes based on tr observations
            gmin1=FSRenvmdr(tr,p,'prob',[0.99; 0.999; 0.01; 0.5],'init',init);
            
            for ii=(i-1):size(gmin1,1)
                
                % CHECK IF STOPPING RULE IS FULFILLED
                % ii>=size(gmin1,1)-2 = final, penultimate or antepenultimate value
                % of the resuperimposed envelope based on tr observations
                if mdr(ii,2)>gmin1(ii,c99) && ii>=size(gmin1,1)-2
                    % Condition S1
                    mes=['$r_{min}('   int2str(mdr(ii,1)) ',' int2str(tr) ')>99$\% envelope'];
                    if msg
                        disp(['Superimposition stopped because r_{min}(' int2str(mdr(ii,1)) ',' int2str(tr) ')>99% envelope']);
                        disp(mes);
                    end
                    sto=1;
                    break
                elseif ii<size(gmin1,1)-2 &&  mdr(ii,2)>gmin1(ii,c999)
                    % Condition S2
                    mes=['$r_{min}('   int2str(mdr(ii,1)) ',' int2str(tr) ')>99.9$\% envelope'];
                    if msg
                        disp(['Superimposition stopped because r_{min}(' int2str(mdr(ii,1)) ',' int2str(tr) ')>99.9% envelope']);
                    end
                    sto=1;
                    break
                else
                    % mdr is inside the envelopes, so keep resuperimposing
                end
            end
            
            if plo==2
                % Plotting part
                subplot(nr,nc,jwind,'Parent',figure2);
                ylim([yl1 yl2]); xlim([xl1 xl2]);
                box('on'); hold('on');
                
                % Show curve of mdr up to step tr-1 (notice that the envelope is
                % based on tr observations. Step tr-1 in matrix mdr is
                % (tr-1)-mdr(1,1)+1=tr-mdr(1,1)
                plot(mdr(1:(tr-mdr(1,1)),1),mdr(1:(tr-mdr(1,1)),2));
                
                % Display the lines associated with 1%, 99% and 99.9% envelopes
                line(gmin1(:,1),gmin1(:,[2 3 4]),'LineWidth',2,'LineStyle','--','Color',[0 0 1]);
                line(gmin1(:,1),gmin1(:,5),'LineWidth',2,'LineStyle','--','Color',[0.3 0.3 0.2]);
                
                strtemp=['$r_{min}(m,' int2str(tr) ')$'];
                % get the position of the current pane
                gposcurax=get(gca,'position');
                % set the width and the height of the current pane but keep
                % unaltered the distance from bottom left corner
                set(gca,'position',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3)*1.3,gposcurax(4)*1.3]);
                % Subplots located on the left hand side
                plleft=1:nc:(nr*(nc-1)+1);
                % Subplots located on the right hand side
                plright=nc:nc:nr*nc;
                % For all plots not located on the left and the right hand side
                % delete numbers on the y axis (YTickLabels)
                getYTickLab=get(gca,'YTickLabel');
                if isempty(intersect(jwind,[plleft plright])) && nr*nc>1
                    set(gca,'YTickLabel',[]);
                end
                % For all plots not located on the right hand side put the
                % yaxis location on the right
                if ~isempty(intersect(jwind,plright)) && nr*nc>1  && nc>1
                    set(gca,'YAxisLocation','right');
                end
                % For all plots which are not on the bottom hand side delete
                % numbers on the x axis (XTickLabels)
                if isempty(intersect(jwind,(nr-1)*nc+1:(nr*nc)))
                    set(gca,'XTickLabel',[]);
                else
                    % For the plots on the bottom side add the xlabel
                    xlabel('Subset size m');
                end
                annotation(figure2,...
                    'textbox',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3),gposcurax(4)+0.05],...
                    'String',{strtemp},'Tag','mes1',...
                    PrVaCell{:});
                
                hold('off');
                jwind=jwind+1;
            end
            
            if sto==1
                if plo==2
                    % Write on the plot the reason why the procedure stopped
                    annotation(figure2,...
                        'textbox',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3),gposcurax(4)],...
                        'String',{mes},'Tag','mes2',...
                        PrVaCell{:});
                    
                    % Unless for all plots not located on the right hand side
                    % For the final plot put the yaxis location on the right
                    % Unless it is the first on the left hand side
                    if isempty(intersect(jwind-1,plleft)) && nr*nc>1  && nc>1
                        set(gca,'YTickLabel',getYTickLab);
                        set(gca,'YAxisLocation','right');
                    end
                    
                    % add X label again to the last plot
                    getaxes=get(figure1,'Children');
                    getaxes=getaxes(end);
                    set(gca,'XTickLabel',get(getaxes,'XTick'))
                    xlabel('Subset size m');
                    
                    % if jwind=2 than the width of the last plot is enlarged to
                    % full screen
                    if jwind==2 && nr*nc>1
                        set(gca,'Position',[0.1 0.1 0.85 0.85])
                        % Relocate the two messages
                        set(findall(gcf,'Tag','mes1'),'Units','normalized','Position',[0.3 0.7 0.5 0.1])
                        set(findall(gcf,'Tag','mes2'),'Units','normalized','Position',[0.3 0.6 0.5 0.1])
                    end
                end
                
                break
            end
            if plo==2
                if jwind==nr*nc+1
                    jwind=1;
                    figure2 = figure('PaperSize',[20.98 29.68],'Name',['Resuperimposed envelopes #' int2str(resup)]);
                    resup=resup+1;
                    % Create axes (Following line should not be necessary
                    %axes('Parent',figure2);
                end
            end
            
        end
        
        %% Stage 2a: subset validation
        % In this part we check whether the subset is homogeneous. In other
        % words we verify conditions H1 and H2
        % tr= m^\dagger+k+1
        % m^\dagger+k=tr-1
        % m*=mdr(ii,1)
        % Condition H2
        % Check if stopping rule takes place at m* <m^\dagger+k
        if (mdr(ii,1)<tr-1)
            % Condition H2b and H2a
            if sum(gmin1(ii+1:end,4)>mdr(ii+1:size(gmin1,1),2))>0
                if msg
                    disp(['Subsample of ' int2str(tr-1) ' units is not homogeneous because the curve was above 99.99% and later it was below 1%']);
                    disp('----------------------------------------');
                end
                % Find m^{1%} that is the step where mdr goes below the 1%
                % threshold for the first time
                % gfind = concatenate all the steps from m^*+1 to m^\dagger+k-1
                gfind=[gmin1(ii+1:end,1) gmin1(ii+1:end,4)>mdr(ii+1:size(gmin1,1),2)];
                % select from gfind the steps where mdr was below 1% threshold
                % gfind(1,1) contains the first step where mdr was below 1%
                gfind=gfind(gfind(:,2)>0,1);
                % find maximum in the interval m^\dagger=mdr(i,1) to the step
                % prior to the one where mdr goes below 1% envelope
                if length(gfind)==1
                    tr=gfind;
                else
                    % Make sure that gfind(1,1)-mdr(1,1) is not smaller
                    % than i otherwise tr  becomes empty
                    if gfind(1,1)-mdr(1,1)>=i
                        tr=sortrows(mdr(i:gfind(1,1)-mdr(1,1),1:2),2);
                        tr=tr(end,1);
                    else
                        tr=mdr(i,1);
                    end
                end
                if msg
                    disp('Probably there are two overlapping groups');
                    disp(['Using the criterion of the maximum, the group of homogenous obs. is=' int2str(tr)]);
                end
                % tr is redefined and is associated with the step associated to
                % the maximum value of r_min
                % try=sormcl[rows(sormcl),1]+1;
                tr=tr+1;
            else
                if msg
                    disp(['Subsample of ' int2str(tr-1) ' units is homogeneous']);
                end
            end
        else
        end
        ndecl=n-tr+1;
        
    elseif  isempty(bonflev) && signal==2 % exceedance of the lower threshold
        
        
        if plo==2
            nr=3; nc=3;
            resup=1;
            jwind=1;
            figure2 = figure('PaperSize',[20.98 29.68],'Name',['Resuperimposed envelopes #' int2str(resup)]);
            
        end
        
        
        
        for tr=(mdag-1):-1:max(mdag-80,init+1)
            % Compute theoretical envelopes based on tr observations
            gmin1=FSRenvmdr(tr,p,'prob',[0.99; 0.01; 0.5; 0.001],'init',init);
            
            ii=size(gmin1,1);
            
            % CHECK IF STOPPING RULE IS FULFILLED
            % resuperimposed envelope based on tr observations is greater
            % than 0.1% threshold
            if mdr(ii,2)>gmin1(ii,5)
                % Condition N1
                mes=['$r_{min}('   int2str(mdr(ii,1)) ',' int2str(tr) ')>0.1$\% envelope'];
                if msg
                    disp(['Superimposition stopped because r_{min}(' int2str(mdr(ii,1)) ',' int2str(tr) ')>1% envelope']);
                    disp(mes);
                end
                sto=1;
                
            elseif tr==max(mdag-80,init+1)
                disp('Exceedance of the lower envelope in the intial part of the search')
                disp('Please decrease the value of init')
                % mdr is still below the 0.1% lower envelope, so keep resuperimposing
            else
            end
            
            
            if plo==2
                % Plotting part
                subplot(nr,nc,jwind,'Parent',figure2);
                ylim([yl1 yl2]); xlim([xl1 xl2]);
                box('on'); hold('on');
                
                % Show curve of mdr up to step tr-1 (notice that the envelope is
                % based on tr observations. Step tr-1 in matrix mdr is
                % (tr-1)-mdr(1,1)+1=tr-mdr(1,1)
                plot(mdr(1:(tr-mdr(1,1)),1),mdr(1:(tr-mdr(1,1)),2));
                
                % Display the lines associated with 1%, 50% and 99% envelopes
                line(gmin1(:,1),gmin1(:,2:4),'LineWidth',2,'LineStyle','--','Color',[0 0 1]);
                
                % Display the lines associated with 0.1% envelopes
                line(gmin1(:,1),gmin1(:,5),'LineWidth',2,'LineStyle','--','Color',[1 0 0]);
                
                
                strtemp=['$r_{min}(m,' int2str(tr) ')$'];
                % get the position of the current pane
                gposcurax=get(gca,'position');
                % set the width and the height of the current pane but keep
                % unaltered the distance from bottom left corner
                set(gca,'position',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3)*1.3,gposcurax(4)*1.3]);
                % Subplots located on the left hand side
                plleft=1:nc:(nr*(nc-1)+1);
                % Subplots located on the right hand side
                plright=nc:nc:nr*nc;
                % For all plots not located on the left and the right hand side
                % delete numbers on the y axis (YTickLabels)
                getYTickLab=get(gca,'YTickLabel');
                if isempty(intersect(jwind,[plleft plright])) && nr*nc>1
                    set(gca,'YTickLabel',[]);
                end
                % For all plots not located on the right hand side put the
                % yaxis location on the right
                if ~isempty(intersect(jwind,plright)) && nr*nc>1  && nc>1
                    set(gca,'YAxisLocation','right');
                end
                % For all plots which are not on the bottom hand side delete
                % numbers on the x axis (XTickLabels)
                if isempty(intersect(jwind,(nr-1)*nc+1:(nr*nc)))
                    set(gca,'XTickLabel',[]);
                else
                    % For the plots on the bottom side add the xlabel
                    xlabel('Subset size m');
                end
                annotation(figure2,...
                    'textbox',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3),gposcurax(4)+0.05],...
                    'String',{strtemp},'Tag','mes1',...
                    PrVaCell{:});
                
                hold('off');
                jwind=jwind+1;
            end
            
            if sto==1
                if plo==2
                    % Write on the plot the reason why the procedure stopped
                    annotation(figure2,...
                        'textbox',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3),gposcurax(4)],...
                        'String',{mes},'Tag','mes2',...
                        PrVaCell{:});
                    
                    % Unless for all plots not located on the right hand side
                    % For the final plot put the yaxis location on the right
                    % Unless it is the first on the left hand side
                    if isempty(intersect(jwind-1,plleft)) && nr*nc>1  && nc>1
                        set(gca,'YTickLabel',getYTickLab);
                        set(gca,'YAxisLocation','right');
                    end
                    
                    % add X label again to the last plot
                    set(gca,'XTickLabel',get(get(figure1,'Children'),'XTick'))
                    xlabel('Subset size m');
                    
                    % if jwind=2 than the width of the last plot is enlarged to
                    % full screen
                    if jwind==2 && nr*nc>1
                        set(gca,'Position',[0.1 0.1 0.85 0.85])
                        % Relocate the two messages
                        set(findall(gcf,'Tag','mes1'),'Units','normalized','Position',[0.3 0.7 0.5 0.1])
                        set(findall(gcf,'Tag','mes2'),'Units','normalized','Position',[0.3 0.6 0.5 0.1])
                    end
                end
                
                break
            end
            if plo==2
                if jwind==nr*nc+1
                    jwind=1;
                    resup=resup+1;
                    figure2 = figure('PaperSize',[20.98 29.68],'Name',['Resuperimposed envelopes #' int2str(resup)]);
                    % Create axes (Following line should not be necessary
                    %axes('Parent',figure2);
                end
            end
            
        end
        
        ndecl=n-tr;
    else
        ndecl=n-mdr(i,1);
    end
    
    if msg
        disp('----------------------------');
        disp('Final output');
        disp(['Number of units declared as outliers=' int2str(ndecl)]);
    end
else
    if msg
        disp('Sample seems homogeneous, no outlier has been found');
    end
    ndecl=0;
end


%% End of the forward search
if msg && isempty(bonflev)
    disp('Summary of the exceedances');
    disp(nout);
end

if msg
    if ~isempty(extram3)
        disp(extram3);
    else
    end
    
    if ~isempty(extram2)
        disp(extram2);
    else
    end
end

group=ones(n,1);

% Plot entry order of the untis
% plot([Un(1,1)-1;Un(:,1)],bb','x')

if ndecl>0
    % Now find the list of the units declared as outliers
    % bsel=~isnan(bb(:,tr-init+1));
    % ListOut=setdiff(1:n,bsel,1);
    % REMARK: If the units forming subset have not been stored for all
    % steps of the fwd search then it is necessary to call procedure FSRbsb
    % to find unit forming subset in step n-decl
    if size(bb,2)<n-init+1
        % then it is necessary to understand what are the units belonging to
        % subset in step n-ndecl.
        % colbb is the column number of bb which contains the units forming
        % subset in the largest step which has been stored among those
        % which are smaller or equal than n-decl.
        % The units in column colbb of matrix bb will form the initial
        % subset in the call of routine FSRbsb
        colbb=find(sum(isnan(bb),1)>=ndecl,1,'last');
        
        %  if sum(~isnan(bb(:,colbb)))<n-ndecl then it is necessary to call
        %  procedure FSRbsb or FSRHbsb or FSRBbsb
        if sum(~isnan(bb(:,colbb)))<n-ndecl
            bsb=seq(~isnan(bb(:,colbb)));
            if isempty(model)
                % Call procedure FSRbsb
                [Un,BB] = FSRbsb(y,X,bsb,'intercept',intercept,'init',n-ndecl,'nocheck',1);
            elseif strcmp(model,'H')
                % Call procedure FSRHbsb
                Z=INP.Z;
                gridsearch=options.gridsearch;
                modeltype=options.modeltype;
                [Un,BB] = FSRHbsb(y,X,Z,bsb,'intercept',intercept,'init',n-ndecl,...
                    'gridsearch',gridsearch,'modeltype',modeltype,'nocheck',1,'msg',0);
            elseif strcmp(model,'B')
                % Call procedure FSRBbsb
                beta0=INP.beta0;
                R=INP.R;
                tau0=INP.tau0;
                n0=INP.n0;
                [Un,BB] = FSRBbsb(y, X, beta0, R, tau0, n0,'bsb',bsb,'intercept',intercept,...
                    'init',n-ndecl,'nocheck',1,'msg',0);
            else
                error('FSDA:FSRcore:WrongModel','Specified model is not supported: possible values are ''H'' (heteroskedastic model) , '''' empty value for linear regression,  ''B'' (Bayesian linear regression)')
            end
            % The first column of BB contains the units forming subset in
            % step n-ndecl
            ListOut=setdiff(seq,BB(:,1));
        else
            ListOut=seq(isnan(bb(:,colbb)));
        end
    else
        ListOut=seq(isnan(bb(:,end-ndecl)));
    end
    % Add to ListOut all the units which have equal values in terms of X
    % and to y to those declared as outliers
    add=zeros(round(n*5),1);
    good=setdiff(seq,ListOut);
    Xy=[X y];
    ij=0;
    for i=1:length(ListOut)
        for j=1:length(good)
            if isequal(Xy(good(j),:),Xy(ListOut(i),:))
                ij=ij+1;
                add(ij)=good(j);
            end
            %   disp(['i' num2str(i) 'j' num2str(j)])
        end
    end
    if ij>0
        add=add(1:ij);
        add=unique(add);
        ListOut=[ListOut,add'];
    end
    % Store the values of beta coefficients in step n-ndecl
    ndecl=length(ListOut);
    
    % Remark: if ndecl>n-init then the number of outliers is set to n-init
    ndecl=min(n-init,ndecl);
    
    if strcmp(model,'H')
        Hetero=INP.Hetero;
        hetero=Hetero(end-ndecl,2:end);
    end
    % Note that the S2 which comes out from procedure FSRBmdr (differently
    % from FSRmdr and FSRHmdr is rescaled, that is it has been inflated by
    % the consistency term, therefore in order to find the unrescaled one
    % we have to recal procedure regressB). Given that the estimate of beta
    % is also affected we also need to recompute it
    if strcmp(model,'B')
        goodobs=setdiff(seq,ListOut);
        beta0=INP.beta0;
        R=INP.R;
        tau0=INP.tau0;
        n0=INP.n0;
        if ~isempty(goodobs)
            outregrB = regressB(y(goodobs), X(goodobs,:), beta0, R, tau0, n0,'intercept',intercept,...
                'nocheck',1);
        else
            % posterior values are equal to prior values, if goodobs is
            % empty
            outregrB=struct;
            outregrB.beta1=beta0;
            outregrB.tau1=tau0;
        end
        scale=sqrt(1/outregrB.tau1);
        beta = outregrB.beta1';
    else
        beta = Bcoeff(end-ndecl,2:end);
        scale= sqrt(S2(end-ndecl,2));
    end
    group(ListOut)=2;
else
    % No outlier is found.
    % Store the values of beta coefficients in final step of the fwd search
    beta = Bcoeff(end,2:end);
    scale= sqrt(S2(end,2));
    ListOut=NaN;
    if strcmp(model,'H')
        Hetero=INP.Hetero;
        hetero=Hetero(end,2:end);
    end
end


%% Scatter plot matrix with the outliers shown with a different symbol

if plo==1 || plo==2
    figure;
    if isempty(options.namey)
        namey=char('y');
    else
        namey=options.namey;
    end
    
    if intercept==1
        if isempty(options.nameX)
            nameX=cellstr(num2str((1:p-1)','X%d'));
        else
            nameX=options.nameX;
        end
        [H,AX,BigAx]=gplotmatrix(X(:,2:end),y,group,'br','+o',[],[],[],nameX,namey);
    else
        if isempty(options.nameX)
            nameX=cellstr(num2str((1:p)','X%d'));
        else
            nameX=options.nameX;
        end
        [H,AX,BigAx]=gplotmatrix(X,y,group,'br','+o',[],[],[],nameX,namey);
    end
    set(gcf,'Name','Scatter plot matrix y|X with outliers highlighted');
    
    % The second condition is necessary because in the Bayesian case all
    % units can be declared as outliers
    if ndecl>0 && n-ndecl>0
        set(H(:,:,1),'DisplayName','Good units');
        set(H(:,:,2),'DisplayName','Outliers');
        % save the indices of the outliers (ListOut) to the
        % 'UserData' field of the second group of H(:,:,2)
        set(H(:,:,2), 'UserData' , ListOut');
    elseif ndecl>0
        set(H(:,:,1),'DisplayName','Outliers');
    else
        set(H,'DisplayName','Units');
    end
    
    % save the indices of the outliers (ListOut) to the
    % 'UserData' field of the second group of H(:,:,2)
    %     if ~isnan(ListOut)
    %         set(H(:,:,2), 'UserData' , ListOut');
    %     end
    
    % The following line adds objects to the panels of the yX
    % add2yX(H,AX,BigAx,outadd,group,ListOut,bivarfit,multivarfit,labeladd)
    add2yX(H,AX,BigAx,'intercept',intercept,'bivarfit',bivarfit,'multivarfit',multivarfit,'labeladd',labeladd);
end

%% Structure returned by function FSR
out=struct;
out.ListOut  = ListOut;
out.outliers = ListOut;


% If you wish that the output also contains the list of units not declared
% as outliers, please uncomment the two following lines.
% ListIn=seq(~isnan(bb(:,end-ndecl)));
% out.ListIn=ListIn;

out.mdr=mdr;
out.Un=Un;
if isempty(bonflev)
    out.nout=nout;
end
out.beta=beta';
if strcmp(model,'H')
    out.hetero=hetero;
end

out.scale=scale;

%% Callback functions used to "pin" quantile labels and vertical line to axes.

    function fprecallback(obj,~)
        % When the panning action starts, the quantile labels and the
        % vertical line which identifies the final part of the search are
        % made invisible.
        
        hanno = findall(obj,'Tag', 'quantile_label');
        set(hanno,'Visible','off');
        hanno2 = findall(obj, 'Tag' , 'FinalPartLine');
        set(hanno2,'Visible','off');
        
    end

    function fpostcallback(obj,evd)
        %When the panning action terminates, the new positions of the
        %labels and of the vertical line are set.
        
        % axis limits
        xlimits = get(evd.Axes,'XLim'); xmin = xlimits(1); xmax=xlimits(2);
        
        % QUANTILES ANNOTATION: the handles
        hanno1 = findall(obj, 'Tag' , 'quantile_label');
        % the quantiles values at each step of the search
        xy = get(hanno1,'UserData');
        % the step 'init'
        init = xy{1,1}(1);
        % the max between init and the lower axis limit
        x = floor(max(xmin , init));
        % the position in the quantiles matrix of x
        xi = x - init + 1;
        % xp is the actual position of the annotation in the plot: normally
        % xp = x, but we move it towards right a bit (by 3 steps) when x is
        % so that the annotation would overlap with the quantiles values
        % labels.
        if x-xmin>0
            xp = x;
        else
            xp = x-3;
        end
        % the figure coordinates which correspond to the anotations positions.
        [y99999, y9999 , y999 , y50 , y99 , y1] = ...
            deal(xy{1,1}(xi,2) , xy{2,1}(xi,2) , xy{3,1}(xi,2) , xy{4,1}(xi,2) , xy{5,1}(xi,2) , xy{6,1}(xi,2));
        [figx, figy1] = dsxy2figxy(evd.Axes, xp, y1);
        [figx, figy99] = dsxy2figxy(evd.Axes, xp, y99);
        [figx, figy50] = dsxy2figxy(evd.Axes, xp, y50);
        [figx, figy999] = dsxy2figxy(evd.Axes, xp, y999);
        [figx, figy9999] = dsxy2figxy(evd.Axes, xp, y9999);
        [figx, figy99999] = dsxy2figxy(evd.Axes, xp, y99999);
        positionValCell(1,1) = {[figx figy99999 0 0]};
        positionValCell(2,1) = {[figx figy9999 0 0]};
        positionValCell(3,1) = {[figx figy999 0 0]};
        positionValCell(4,1) = {[figx figy50 0 0]};
        positionValCell(5,1) = {[figx figy99 0 0]};
        positionValCell(6,1) = {[figx figy1 0 0]};
        % set the new figure coordinates and make the annotations visible again
        set(hanno1 , {'Position'} , positionValCell, 'Visible','on');
        
        % Uncomment the next lines to display in the command window the new
        % positions
        % disp(['Min axis value:' int2str(xmin)]);
        % disp(['x value:' int2str(x)]);
        % disp(['xi value:' int2str(xi)]);
        % disp(['xp value:' int2str(xp)]);
        
        % VERTICAL LINE which identifies the final part of the search.
        % For that line we only need to recompute the X position.
        hanno2 = findall(obj, 'Tag' , 'FinalPartLine');
        LineData = get(hanno2,'UserData');
        istep = LineData(1); yl1 = LineData(2); yl2 = LineData(3);
        if istep > xmax, istep = xmax; end
        if istep < xmin, istep = xmin; end
        [figx, figy]  = dsxy2figxy(evd.Axes, istep, yl1);
        
        positionLineOri=get(hanno2,'Position');
        positionLineOri(1)=figx;
        
        positionLineCell = {positionLineOri};
        xCell = {[figx figx]};
        
        set(hanno2 , {'Position'} , positionLineCell, {'X'} , xCell, 'Visible','on');
        
    end
end
