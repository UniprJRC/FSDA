function out=RobCov(X,scaledres,scaleest,varargin)
%robCov computes covariance matrix of robust regression coefficients
%
%<a href="matlab: docsearchFS('RobCov')">Link to the help function</a>
%
% Under some regularity conditions, robust (S and MM) estimates are
% asymptotically normal, thereby allowing for Wald-type tests and
% confidence intervals. The covariance matrix of the estimated parameters
% \[
%   cov(\hat \beta)= q^2 \times \sigma^2 \times v \times V_X^{-1}
% \]
% consists of four parts:
%1) $q$ a correction factor for the scale estimate;
% 2) $\sigma$ the scale parameter.
% 3) $v$ a correction factor depending on the $\psi$ function which is
% used;
% 4) $V_X$= a matrix part. For OLS $V_X=X'X$. Given that in robust
% regression we give a weight to each observation, the matrix $X'X$ should
% be replaced by something like $X'WX$, where $W$ is a diagonal matrix
% containing the weights assigned to each observation.
% The purpose of this function is to provide the user with different
% options for the estimate of $cov(\hat \beta)$ where $\hat \beta$ is  a
% vector of regression coefficients obtained using S or MM estimation and a
% particular $\rho$ function.
%
%
%
%  Required input arguments:
%
%    X :     Data matrix of explanatory variables (also called 'regressors')
%            of dimension (n x p-1). Rows of X represent observations, and
%            columns represent variables.
%               Data Types - single | double
% scaledres : Scaled residuals.Vector. n-times-1 vector containing scaled
%             residuals $r_i/\hat \sigma$.
%               Data Types - single | double
% scaleest  : robust estimate of the scale. Scalar. Robust estimate of
%             sigma ($\hat \sigma$).
%               Data Types - single | double
%
%  Optional input arguments:
%
%  intercept :  Indicator for constant term. Scalar. If 1, a model with
%               constant term will be fitted (default), else no constant
%               term will be included.
%               Example - 'intercept',1
%               Data Types - double
%      eff     : nominal efficiency. Scalar.
%                Scalar defining nominal efficiency (i.e. a number between
%                 0.5 and 0.99). The default value is 0.95
%                 Asymptotic nominal efficiency is:
%                 $(\int \psi' d\Phi)^2 / (\psi^2 d\Phi)$
%                 Example - 'eff',0.99
%                 Data Types - double
%  intercept :  Indicator for constant term. Scalar. If 1, a model with
%               constant term will be fitted (default), else no constant
%               term will be included.
%               Example - 'intercept',1
%               Data Types - double
%         bdp :  breakdown point. Scalar.
%               It measures the fraction of outliers
%               the algorithm should resist. In this case any value greater
%               than 0 but smaller or equal than 0.5 will do fine.
%               Note that given bdp nominal
%               efficiency is automatically determined.
%                 Example - 'bdp',0.4
%                 Data Types - double
%               REMARK: just one between bdp and eff must be specified. If
%               both of them are specified an error is produced. If both of
%               them are not specified the defulat is tu use the tuning
%               constant associated to a nominal efficiency of 0.95.
%     rhofunc : rho function. String. String which specifies the rho function which must be used to
%               weight the residuals. Possible values are
%               'bisquare'
%               'optimal'
%               'hyperbolic'
%               'hampel'.
%               'bisquare' uses Tukey's $\rho$ and $\psi$ functions.
%               See TBrho.m and TBpsi.m.
%               'optimal' uses optimal $\rho$ and $\psi$ functions.
%               See OPTrho.m and OPTpsi.m.
%               'hyperbolic' uses hyperbolic $\rho$ and $\psi$ functions.
%               See HYPrho.m and HYPpsi.m.
%               'hampel' uses Hampel $\rho$ and $\psi$ functions.
%               See HArho.m and HApsi.m.
%               The default is bisquare
%                 Example - 'rhofunc','optimal'
%                 Data Types - double
% rhofuncparam: Additional parameters for the specified rho function.
%               Scalar or vector.
%               For hyperbolic rho function it is possible to set up the
%               value of k = sup CVC (the default value of k is 4.5).
%               For Hampel rho function it is possible to define parameters
%               a, b and c (the default values are a=2, b=4, c=8)
%                 Example - 'rhofuncparam',5
%                 Data Types - single | double
%  Output:
%
%  out :     A structure containing the following fields
%
%   out.covrob = p-times-p (if intercept is 1 else is (p-1)-by-(p-1)) matrix
%               containing asymptotic variance covariance
%               matrix of regression coefficients. covrob implements
%               equation (4.49) of p. 101 of Maronna et al. (2006)
%               namely:
%                \[
%                \mbox{covrob} = cov( \hat \beta) = \hat \sigma^2 \hat v (X'X)^{-1}
%                \]
%                where
%                \[
%                \hat v =  \frac{n}{n-p} n\frac{\sum_{i=1}^n \psi(r_i/\hat \sigma)^2}{\sum_{i=1}^n \psi'(r_i/\hat \sigma)^2}
%                \]
%  out.covrob1 =  p-times-p (if intercept is 1 else is (p-1)-by-(p-1)) matrix
%               containing asymptotic variance covariance
%               matrix of regression coefficients. covrob1 implements
%               equation (7.81) of p. 171 of Huber and Ronchetti (2009)
%               with $(X'X)^{-1}$ replaced by $(X' W X)^{-1}$
%               namely:
%                \[
%                 \mbox{covrob1} =  K^2  \hat v  (X' W X)^{-1};
%                \]
%                 where $K=1+p n \frac{var(\psi' (r/\hat \sigma))}{ \left[ (\sum_{i=1}^n
%                 \psi'(r_i/\hat \sigma)) \right]^2}$ ;
%  out.covrob2 =  p-times-p (if intercept is 1 else is (p-1)-by-(p-1)) matrix
%               containing asymptotic variance covariance
%               matrix of regression coefficients. covrob1 implements
%               equation (7.81) of p. 171 of Huber and Ronchetti (2009)
%               with $X'X$ and $K^2$
%               namely:
%                \[
%                 \mbox{covrob2} =  K^2  \hat v  (X' X)^{-1};
%                \]
%                 where $K=1+p n \frac{var(\psi' (r/\hat \sigma))}{ (\sum_{i=1}^n
%                 \psi'(r_i/\hat \sigma))^2}$ ;
%  out.covrob3 =  p-times-p (if intercept is 1 else is (p-1)-by-(p-1)) matrix
%               containing asymptotic variance covariance
%               matrix of regression coefficients. covrob implements
%               equation (7.82) of p. 171 of of Huber and Ronchetti (2009).
%               namely:
%                \[
%                \mbox{covrob3} =  K  \hat v  (X' W X)^{-1};
%                \]
%                 where $K=1+p n \frac{var(\psi' (r/\hat \sigma))}{ (\sum_{i=1}^n
%                 \psi'(r_i/\hat \sigma))^2}$ ;
%  out.covrob4 =  p-times-p (if intercept is 1 else is (p-1)-by-(p-1)) matrix
%               containing asymptotic variance covariance
%               matrix of regression coefficients. covrob implements
%               equation (7.83) of p. 171 of of Huber and Ronchetti (2009).
%               namely:
%                \[
%                 \mbox{covrob4} =  \frac{1}{n-p} K^{-1} \sum_{i=1}^n (\psi(r_i/\hat \sigma))^2  (X' W X)^{-1} X'X (X' W X)^{-1};
%                \]
%      out.q =  scalar. Correction for scale estimate (see Maronna and Yohai
%               CSDA 2010). It is defined as
%               \[
%               q=1+\frac{p}{2n} \frac{a}{b \times c}
%               \]
%               where
%               \[
%               a = \frac{1}{n} \sum_{i=1}^n (\psi(r_i/\hat \sigma))^2
%               \]
%               \[
%               b = \frac{1}{n} \sum_{i=1}^n (\psi'(r_i/\hat \sigma))^2
%               \]
%               \[
%               c = \frac{1}{n} \sum_{i=1}^n (\psi(r_i/\hat \sigma)) r_i/\hat \sigma)
%               \]
%
% See also: Sreg, MMreg, Taureg
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
% Huber, P.J. and Ronchetti, E.M. (2009), "Robust Statistics, 2nd Edition",
% Wiley.
% Maronna, R.A., and Yohai V.J. (2010), Correcting MM estimates for fat data
% sets, "Computational Statistics and Data Analysis", Vol. 54, pp. 3168-3173.
% Koller, M. and W. A. Stahel (2011), Sharpening wald-type inference in
% robust regression for small samples, "Computational Statistics & Data
% Analysis", Vol. 55, pp. 2504-2515.
% Croux, C., G. Dhaene, and D. Hoorelbeke (2003). Robust standard errors
% for robust estimators. Technical report, Dept. of Applied Economics, K.U.
% Leuven.
%
% Copyright 2008-2018.
% FSDA toolbox
%
%
%<a href="matlab: docsearchFS('RobCov')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Compare the 5 estimates of cov matrix.
    n=200;
    p=3;
    state1=123456;
    randn('state', state1);
    X=randn(n,p);
    y=randn(n,1);
    kk=10;
    ycont = y;
    ycont(1:kk)=ycont(1:kk)+7;
    [outS]=Sreg(ycont,X);
    rhofunc='optimal';
    bdp=0.5;
    out=RobCov(X,outS.residuals,outS.scale);
    disp('Compare 5 estimates of cov(\hat beta)')
    disp(out.covrob)
    disp('--------')
    disp(out.covrob1)
    disp('--------')
    disp(out.covrob2)
    disp('--------')
    disp(out.covrob3)
    disp('--------')
    disp(out.covrob4)
%}

%{
    % Compare t stat from S and MM estimator.
    rhofunc='optimal';
    bdp=0.5;
    out=RobCov(X,outS.residuals,outS.scale,'rhofunc',rhofunc,'bdp',0.5);
    covrobS=out.covrob;
    covrobS1=out.covrob1;
    covrobS2=out.covrob2;
    covrobS3=out.covrob3;
    covrobS4=out.covrob4;

    % Compute robust S t-statistics
    tstatS=outS.beta./sqrt(diag(covrobS));
    tstatS1=outS.beta./sqrt(diag(covrobS1));
    tstatS2=outS.beta./sqrt(diag(covrobS2));
    tstatS3=outS.beta./sqrt(diag(covrobS3));
    tstatS4=outS.beta./sqrt(diag(covrobS4));

    eff=0.95;
    outMM=MMregcore(ycont,X,outS.beta,outS.scale);
    out=RobCov(X,outMM.residuals,outS.scale,'rhofunc',rhofunc,'eff',eff);
    covrobMM=out.covrob;
    covrobMM1=out.covrob1;
    covrobMM2=out.covrob2;
    covrobMM3=out.covrob3;
    covrobMM4=out.covrob4;
    tstatMM=outMM.beta./sqrt(diag(covrobMM));
    tstatMM1=outMM.beta./sqrt(diag(covrobMM1));
    tstatMM2=outMM.beta./sqrt(diag(covrobMM2));
    tstatMM3=outMM.beta./sqrt(diag(covrobMM3));
    tstatMM4=outMM.beta./sqrt(diag(covrobMM4));
    disp('tstat from S')
    disp([tstatS tstatS1 tstatS2 tstatS3 tstatS4])
    disp('--------')
    disp('tstat from MM')
    disp([tstatMM tstatMM1 tstatMM2 tstatMM3 tstatMM4])
    qhat=out.q;
    disp('tstat from MM after correction for sigma')
    disp([tstatMM/qhat tstatMM1/qhat tstatMM2/qhat tstatMM3/qhat tstatMM4/qhat])
%}

%% Beginning of code

% rho (psi) function which has to be used to weight the residuals
rhofuncdef='bisquare';
%rhofuncdef='optimal';

bdpdef='';
effdef='';
% store default values in the structure options
options=struct('intercept',1,'rhofunc',rhofuncdef,'rhofuncparam','','bdp',bdpdef,'eff',effdef);

% check user options and update structure options
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:RobCov:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 3
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

[n,p]=size(X);

intercept=options.intercept;
if intercept
    X = [ones(n,1),X];
    p=p+1;
end


bdp = options.bdp;              % break down point
eff = options.eff;              % nominal efficiency
% Remark: one of them must be empty

if  ~isempty(bdp) && ~isempty(eff)
    error('FSDA:RobCov:WrongInputOpt','Both bdp and eff cannot be specified')
end

if  isempty(bdp) && isempty(eff)
    eff=0.95;
end

rhofunc=options.rhofunc;        % String which specifies the function to use to weight the residuals


% In this case just bdp has been specified
if  ~isempty(bdp)
    if strcmp(rhofunc,'bisquare')
        c=TBbdp(bdp,1);
        psifunc='TB';
        
    elseif strcmp(rhofunc,'optimal')
        % Optimal rho function is strictly increasing on [0 3c] and constant (equal to 3.25c^2) on [3c \infty)
        % E(\rho) = kc = (3.25c^2)*bdp = TBrho(3*c,c)*bdp, being kc the K of
        % Rousseeuw and Leroy (1987)
        
        % Compute tuning constant associated to the requested breakdown
        % point
        % For bdp =0.5 and optimal rho function c1=0.4046
        % Remark: given that in function OPTbdp rho function is defined in the interval 0---2c/3, 2c/3---3c/3, >3c/3
        % it is necessary to divide the output of OPTbdp by 3
        c=OPTbdp(bdp,1)/3;
        psifunc='OPT';
        
    elseif strcmp(rhofunc,'hyperbolic')
        
        if isempty(options.rhofuncparam)
            kdef=4.5;
        else
            kdef=options.rhofuncparam;
        end
        
        % Use (if possible) precalculated values of c,A,b,d and kc
        if kdef == 4 && bdp==0.5
            c =2.158325031399727;
            A =1.627074124322223e-04;
            B =0.006991738279441;
            d =0.016982948780061;
            
        elseif kdef == 4.5 && bdp==0.5
            c =2.010311082005501;
            A =0.008931591866092;
            B =0.051928487236632;
            d =0.132017481327058;
        elseif kdef == 5 && bdp==0.5
            c =1.900709968805313;
            A =0.023186529890225;
            B =0.083526860351552;
            d =0.221246910095216;
        elseif kdef == 4.5 && bdp==0.25
            c =2.679452645778656;
            A =0.464174145115400;
            B =0.588821276233494;
            d =1.092639541625978;
            
        else
            
            % Compute tuning constant associated to the requested breakdown
            % point
            [c,A,B,d]=HYPbdp(bdp,1,kdef);
        end
        
        
        c=[c,kdef,A,B,d];
        psifunc='HYP';
        
        
    elseif strcmp(rhofunc,'hampel')
        
        if isempty(options.rhofuncparam)
            abc=[2,4,8];
        else
            abc=options.rhofuncparam;
        end
        
        % Compute tuning constant associated to the requested breakdown
        % point
        c=HAbdp(bdp,1,abc);
        
        
        c=[c,abc];
        psifunc='HA';
    else
        
        error('FSDA:RobCov:WrongInputOpt','Specified rho function is not supported: possible values are ''bisquare'' , ''optimal'',  ''hyperbolic'', ''hampel''')
    end
end

% In this case just efficiency has been specified
if  ~isempty(eff)
    if strcmp(rhofunc,'bisquare')
        
        c=TBeff(eff,1);
        psifunc='TB';
        
    elseif strcmp(rhofunc,'optimal')
        
        
        % Compute tuning constant associated to the requested nominal efficiency
        % c2 = consistency factor for a given value of efficiency
        % Remark: given that in function OPTeff rho function is defined in the interval 0---2c/3, 2c/3---3c/3, >3c/3
        % it is necessary to divide the output of OPTeff by 3
        c=OPTeff(eff,1)/3;
        
        psifunc='OPT';
        
    elseif strcmp(rhofunc,'hyperbolic')
        
        if isempty(options.rhofuncparam)
            kdef=4.5;
        else
            kdef=options.rhofuncparam;
        end
        
        
        if kdef == 4 && eff==0.85
            c2 =3.212800979614258;
            A2 =0.570183575755717;
            B2 =0.696172437281084;
            d2 =1.205900263786317;
        elseif kdef == 4.5 && eff==0.85
            c2 =3.032387733459473;
            A2 =0.615717108822885;
            B2 = 0.723435958485131;
            d2 =1.321987605094910;
        elseif kdef == 5 && eff==0.85
            c2 =2.911890029907227;
            A2 =0.650228046997054;
            B2 =0.743433840145084;
            d2 =1.419320821762087;
            
        elseif kdef == 4 && eff==0.95
            c2 =4.331634521484375;
            A2 =0.754327484845243;
            B2 =0.846528826589308;
            d2 =1.480099129676819;
        elseif kdef == 4.5 && eff==0.95
            c2 =3.866390228271484;
            A2 =0.791281464739131;
            B2 =0.867016329355630;
            d2 =1.610621500015260;
        elseif kdef == 5 && eff==0.95
            c2 =3.629499435424805;
            A2 =0.818876452066880;
            B2 =0.882004888111327;
            d2 =1.723768949508668;
            
        else
            
            % Compute tuning constant associated to the requested nominal efficiency
            % c2 = consistency factor for a given value of efficiency
            [c2,A2,B2,d2]=HYPeff(eff,1,kdef);
        end
        
        
        c=[c2,kdef,A2,B2,d2];
        psifunc='HYP';
        
    elseif strcmp(rhofunc,'hampel')
        
        if isempty(options.rhofuncparam)
            abc=[2,4,8];
        else
            abc=options.rhofuncparam;
        end
        
        
        % Compute tuning constant associated to the requested nominal efficiency
        % c2 = consistency factor for a given value of efficiency
        c=HAeff(eff,1,abc);
        
        c=[c,abc];
        psifunc='HA';
    else
        error('FSDA:RobCov:WrongInputOpt','Specified rho function is not supported: possible values are ''bisquare'' , ''optimal'',  ''hyperbolic'', ''hampel''')
        
    end
    
end

% c1=TBeff(0.95,1);
% psi=TBpsi(scaledres,c1);
% dpsi=TBpsider(scaledres,c1);

XXpsi=strcat(psifunc,'psi');
XXpsi=str2func(XXpsi);
psi=feval(XXpsi,scaledres,c);

XXpsider=strcat(psifunc,'psider');
XXpsider=str2func(XXpsider);
psider=feval(XXpsider,scaledres,c);


XX=X'*X;
%% Find covrob
% Epsi2=(psi'*psi)/n=  sumpsi2 /n ;
% EXX=(1/n)*XX;
% Edpsi=(1/n)*sum(psider);
% % cov2s should be the one producing traditional standard errors
% covrob=(1/n)*scaleest^2*Epsi2*inv(EXX)/Edpsi^2;
sumpsi2=psi'*psi;
sumpsider=sum(psider);
% vhat implements equation (4.49) of p. 101 of Maronna et al. (2006)
vhat=(n/(n-p))*n*(scaleest^2)*sumpsi2/(sumpsider^2);
invXX=inv(XX);
% Equation below 4.50 of Maronna et al. (2006)
covrob=vhat*invXX; %#ok<MINV>

%% Find covrob1
% See equation 7.81 of Huber and Ronchetti (2009) with X'X replaced by X'WX
% (and sigma is estimated)

% Ksquare =see equation 7.84 of Huber and Ronchetti (2009) P. 171
Ksquare=(1+(p/n)*var(psider,1)/((sumpsider/n)^2))^2;
% disp('varpsider')
% disp(var(psider))
% disp('M^2 varpsider')
% disp(((sumpsider/n)^2))
% disp('psider')
% % disp(psider)
% plot(abs(sort(psider)))
% disp('---------')

K=sqrt(Ksquare);


XXwei=strcat(psifunc,'wei');
XXwei=str2func(XXwei);
w=feval(XXwei,scaledres,c);
w1=sqrt(w);
% Compute (X' W X) / mean(w)
Xw=bsxfun(@times,X,w1);
XWX=(Xw'*Xw)./mean(w);
% Compute ( (X' W X) / mean(w))^-1
invXWX=inv(XWX);

% See equation 7.81 of Huber and Ronchetti (2009) with X'X replaced by X'WX
covrob1=Ksquare*vhat*invXWX;  %#ok<MINV>

%% Find covrob2
% See equation 7.81 of Huber and Ronchetti (2009) with X'X
covrob2=Ksquare*vhat*invXX;

%% Find covrob3
% See equation 7.82 of Huber and Ronchetti (2009)
covrob3=K*(n/(n-p))*(scaleest^2)*(sumpsi2/sumpsider)*invXWX; %#ok<MINV>

%% Find covrob4
% See equation 7.83 of Huber and Ronchetti (2009)
covrob4=(1/(n-p))*(scaleest^2)*sumpsi2*(invXWX*XX*invXWX)/K; %#ok<MINV>

% Find qhat see Maronna and Yohai CSDA 2010 p. 3170 equation (8)
ahat=sumpsi2/n;
bhat=sumpsider/n;
chat=sum(psi.*scaledres)/n;
q=1+(p/(2*n))*(ahat/(bhat*chat));

out=struct;
out.covrob=covrob;
out.covrob1=covrob1;
out.covrob2=covrob2;
out.covrob3=covrob3;
out.covrob4=covrob4;
out.q=q;
end
%FScategory:REG-Regression