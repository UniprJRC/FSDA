function out = RhoPsiWei(u, v, varargin)
%RhoPsiWei finds rho, psi, psi', w functions given bdp, or eff or tuning constant c
%
%<a href="matlab: docsearchFS('RhoPsiWei')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector.
%               Vector containing residuals or Mahalanobis distances
%               for the n units of the sample. If u is empty fields rho, wei,
%               psi, psix and psider of output structure out will be empty.
%               Data Types - single|double
%
%    v         : number of response variables. Scalar. e.g. in regression
%                v=1. So far v can be set to a value greater than 1 just
%                for TB and OPT functions
%               Data Types - single|double|int32|int64
%
%  Optional input arguments:
%
%      bdp    : breakdown point. Scalar. Scalar defining breakdown point
%               (i.e a number between 0 and 0.5). If none of the optional
%               input arguments c, bdp and eff, are not specified bdp is
%               set to 0.5. Note that if c or eff are specified bdp is
%               automatically found.
%               Example - 'bdp',0.2
%               Data Types - single|double
%
%      c      : tuning constant. Scalar. Scalar defining the tuning constant which
%               determines bdp or eff. c corresponds to alpha for mdpd
%               link. Note that if bdp or eff are specified c is
%               automatically found.
%               Example - 'c',3
%               Data Types - single|double
%
%      eff    : asymptotic efficiency. Scalar. Scalar defining requested asymptotic efficiency
%               of the estimate (in general a number in the range 0.5 0.99.
%               Typical values of eff are 0.85, 0.90 or 0.95.
%               Note that if bdp or c are specified eff is
%               automatically found.
%               Example - 'eff',0.9
%               Data Types - single|double
%
%    rhofunc : link function. String. String which specifies the rho
%               function which must be used.
%               Possible values are
%               'TB' or ''biweight', 'OPT' or 'optimal', 'HYP' or 'hyperbolic';
%               'HA' or 'hampel', 'PD' or 'mdpd', 'HU' or 'huber'.
%               'biweight' uses Tukey's $\rho$ and $\psi$ functions.
%               See TBrho and TBpsi.
%               'optimal' uses optimal $\rho$ and $\psi$ functions.
%               See OPTrho and OPTpsi.
%               'hyperbolic' uses hyperbolic $\rho$ and $\psi$ functions.
%               See HYPrho and HYPpsi.
%               'hampel' uses Hampel $\rho$ and $\psi$ functions.
%               See HArho and HApsi.
%               'mdpd' uses Minimum Density Power Divergence $\rho$ and $\psi$ functions.
%               See PDrho.m and PDpsi.m.
%               'hu' uses Hampel $\rho$ and $\psi$ functions.
%               See HArho and HApsi.
%               The default is bisquare
%                 Example - 'rhofunc','optimal'
%                 Data Types - character
%
% rhofuncparam: Additional parameters for the specified rho function.
%               Scalar or vector.
%               For hyperbolic rho function it is possible to set up the
%               value of k = sup CVC (the default value of k is 4.5).
%               For Hampel rho function it is possible to define parameters
%               a, b and c (the default values are a=2, b=4, c=8). For
%               Hyperbolic link function it is possible to pass either one
%               parameter, k (default is k=4.5) or four parameters: k, A, B
%               and d.
%                 Example - 'rhofuncparam',5
%                 Data Types - single | double
%
% Output:
%
%         out:   structure which contains the following fields
%
% out.rho   = vector of length n which contains the rho function associated
%             to the residuals or Mahalanobis distances for the n units of
%             the sample. This field is empty if input u is empty.
% out.psi   = vector of length n which contains the psi function associated
%             to the residuals or Mahalanobis distances for the n units of
%             the sample. This field is empty if input u is empty.
% out.psider = vector of length n which contains the derivative of the
%             psi function associated to the residuals or Mahalanobis
%             distances for the n units of the sample. This field is empty
%             if input u is empty.
% out.wei = vector of length n which contains the weights associated to the
%             residuals or Mahalanobis distances for the n units of the
%             sample. This field is empty if input u is empty.
% out.psix   = vector of length n which contains the psi function mutiplied
%             by u associated to the residuals or Mahalanobis distances for
%             the n units of the sample. This field is empty if input u is
%             empty.
% out.class = Character which specifies the link
%               function which has be used.
%               Possible values are 'bisquare', 'optimal', 'hyperbolic';
%               'hampel', 'huber' or 'mdpd'.
% out.bdp     = scalar which contains the bdp which has been used.
% out.eff     = scalar which contains the eff which has been used.
% out.c1      =  consistency factor (and other parameters)
%                   associated to required breakdown point or nominal
%                   efficiency.
%                   More precisely, out.c1(1) contains consistency
%                   factor c associated to required breakdown point or
%                   nominal efficiency out.c1(2:end) contain other
%                   parameters associated with the rho (psi) function.
%                   For example, if rhofunc='hampel', c1(2:4) must
%                   contain parameters (a, b and c) of Hampel rho function.
%                   If rhofunc is hyperbolic out.c1(1) specifies the value of the tuning
%                constant c (scalar greater than 0 which controls the
%                robustness/efficiency of the estimator)
%                out.c1(2) contains the prefixed value k (sup of the
%                change-of-variance sensitivity) and out.c1(3:5)
%                contain parameters A, B and d.
%   out.kc1   =   Expectation of rho associated with out.c1(1) to get a
%                 consistent estimator at the model distribution
%                 kc1=E(rho)=sup(rho)*bdp
%
%
% See also: TBrho, OPTeff, HYPc, PDeff, HUeff
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('RhoPsiWei')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

%
% Examples:
%
%{
    %% Example of call to RhoPsiWei with 2 input arguments.
    % In this case TB link is invoved with bdp=0.5;
    n=20;
    u=randn(n,1);
    out=RhoPsiWei(u,1);
%}

%{
    %% Example of call to RhoPsiWei with bdp prespecified.
    % In this case TB link is invoved with the bdp specified as input
    n=20;
    u=randn(n,1);
    out=RhoPsiWei(u,1,'bdp',0.1);
%}

%{
    %% Example of call to RhoPsiWei with bdp and rhofunc prespecified.
    % In this case TB link is invoved with the bdp specified as input
    n=20;
    u=randn(n,1);
    rhofunc='PD';
    out=RhoPsiWei(u,1,'bdp',0.1,'rhofunc',rhofunc);
%}

%{
    %% Example of call to RhoPsiWei with eff and rhofunc prespecified.
    % In this case TB link is invoved with the bdp specified as input
    n=20;
    u=randn(n,1);
    rhofunc='OPT';
    eff=0.9;
    out=RhoPsiWei(u,1,'eff',eff,'rhofunc',rhofunc);
%}

%{
    %% Example of call to RhoPsiWei with c prespecified.
    % In this case TB link is invoved with the bdp specified as input
    n=20;
    u=randn(n,1);
    rhofunc='HA';
    c=1.5;
    out=RhoPsiWei(u,1,'c',c,'rhofunc',rhofunc);
%}

%{
    %% Example of call to RhoPsiWei with c prespecified, and personalized input parameters for HA link.
    % In this case TB link is invoved with the bdp specified as input
    n=20;
    u=randn(n,1);
    rhofunc='HA';
    rhofuncparam=[2 3 8];
    eff=.95;
    c=0.5;
    out=RhoPsiWei(u,1,'eff',eff,'rhofunc',rhofunc,'rhofuncparam',rhofuncparam);
%}

%{
    %% Example of call to RhoPsiWei with eff prespecified, and personalized input parameters fpr HYP link.
    % In this case TB link is invoved with the bdp specified as input
    n=20;
    u=-4:0.01:4;
    rhofunc='HYP';
    k=4.2;
    eff=0.9;
    out=RhoPsiWei(u,1,'eff',eff,'rhofunc',rhofunc,'rhofuncparam',k);
    subplot(2,3,1)
    plot(u,out.rho)
    title('HYP: rho function')
    
    subplot(2,3,2)
    plot(u,out.psi)
    title('HYP: psi function')
    
    subplot(2,3,3)
    plot(u,out.psider)
    title('HYP: derivative of  psi function')
    
    subplot(2,3,4)
    plot(u,out.psider)
    title('HYP: psi(x)*x function')
    
    subplot(2,3,5)
    plot(u,out.wei)
    title('HYP: weight function')
%}

%{
    %% Example of calling RhoPsiWei with all parameters related to the constant: c, k, A, B and d
    % First call RhoPsiWei with efficiency fixed and then use the
    %   calculated constant for a second call - it will be much faster
    n=20;
    u=-4:0.01:4;
    rhofunc='HYP';
    k=4.2;
    eff=0.9;
    out=RhoPsiWei(u,1,'eff',eff,'rhofunc',rhofunc,'rhofuncparam',k);
    out1=RhoPsiWei(u,1,'c',out.c1(1),'rhofunc',rhofunc,'rhofuncparam',out.c1(2:5));
    isequal(out.rho, out1.rho)
    isequal(out.psi, out1.psi)
    [bdp1, eff1]=HYPc(out.c1(1), 1, 'k', out.c1(2), 'param', out.c1(3:5));
    isequal(out.bdp, bdp1)
    isequal(out.eff, eff1)

%}
%{
    % Example of call to RhoPsiWei with both bdp and eff.
    % In this case, given that it is possible to supply just one between bdp,
    % eff and c the call throws an error
    n=20;
    u=-4:0.01:4;
    rhofunc='HYP';
    k=4.2;
    eff=0.85;
    out=RhoPsiWei(u,1,'bdp',0.2,'eff',eff,'rhofunc',rhofunc,'rhofuncparam',k);
%}

%{
    % Example of call to RhoPsiWei with u empty
    % In this case fiels psi, psider, psix and wei are empty
    n=20;
    u=randn(n,1);
    rhofunc='PD';
    eff=0.9;
    out=RhoPsiWei([],1,'eff',eff,'rhofunc',rhofunc);
%}

%% Beginning of code

if nargin > 2
    % Check that maximum one input between bdp, eff and c has been supplied
    UserOptions=varargin(1:2:length(varargin));

    checkbdp = strcmp(UserOptions,'bdp');
    checkeff = strcmp(UserOptions,'eff');
    checkc = strcmp(UserOptions,'c');

    ckinput=any(checkbdp)+any(checkeff)+any(checkc);
    if ckinput==0 % in this case nither bpd nor eff nor c has been specified
        bdpdef=0.5;
    else
        bdpdef=[];
    end

    if ckinput>1
        error('FSDA:RhoPsiWei:WrongInputOpt','Only a value among bdp, eff and c must be supplied');
    end


    options=struct('bdp',bdpdef,'eff','','c',[],'rhofunc','TB','rhofuncparam',[]);

    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:RhoPsiWei:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end


    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end


    bdp=options.bdp;
    eff=options.eff;
    c=options.c;
    rhofunc=options.rhofunc;
    rhofuncparam=options.rhofuncparam;
else
    rhofunc='TB';
    bdp=0.5;
    eff=[];
    c=[];
    rhofuncparam=[];
end

out=struct;
if  isempty(u)
    unonempty=false;
else
    unonempty=true;
end

if strcmp(rhofunc,'TB') || strcmp(rhofunc,'biweight') || strcmp(rhofunc,'bisquare')
    % TUKEY BIWEIGHT
    if ~isempty(c)
        bdp=TBc(c,v);
    elseif ~isempty(eff)
        ceff=TBeff(eff,v);
        bdp=TBc(ceff,v);
    else
    end

    out.class='TB';
    c=TBbdp(bdp,v);
    % kc = E(rho) = sup(rho)*bdp
    kc=c^2/6*bdp;
    out.c1=c;
    out.kc1=kc;
    out.bdp=bdp;

    % Find efficiency
    [~,eff]=TBc(c,v);
    out.eff=eff;

    if unonempty==true
        out.rho=TBrho(u,c);
        out.psi=TBpsi(u,c);
        out.psider=TBpsider(u,c);
        out.psix=TBpsix(u,c);
        out.wei=TBwei(u,c);
    end

elseif strcmp(rhofunc,'OPT') || strcmp(rhofunc,'optimal')
    % OPTIMAL
    if ~isempty(c)
        bdp=OPTc(c,v);
    elseif ~isempty(eff)
        ceff=OPTeff(eff,v);
        bdp=OPTc(ceff,v);
    else
    end

    out.class='OPT';
    c=OPTbdp(bdp,v);
    out.c1=c;
    out.kc1=bdp;
    out.bdp=bdp;

    % Find efficiency
    [~,eff]=OPTc(c,v);
    out.eff=eff;

    if unonempty==true
        out.rho=OPTrho(u,c);
        out.psi=OPTpsi(u,c);
        out.psider=OPTpsider(u,c);
        out.psix=OPTpsix(u,c);
        out.wei=OPTwei(u,c);
    end

elseif strcmp(rhofunc,'HA') || strcmp(rhofunc,'hampel')
    % HAMPEL
    out.class='HA';
    if isempty(rhofuncparam)
        abc=[2 4 8];
    else
        % make sure that abc is a row vector
        abc=rhofuncparam(:)';
    end

    if ~isempty(c)
        bdp=HAc(c,v,'param',abc);
    elseif ~isempty(eff)
        ceff=HAeff(eff,v,abc);
        bdp=HAc(ceff,v,'param',abc);
    else
    end


    % Compute tuning constant associated to the requested breakdown
    % point
    cHA=HAbdp(bdp,v,abc);
    % kc = E(rho) = sup(rho)*bdp
    c1=[cHA abc];
    out.c1=c1;
    % kc = E(rho) = sup(rho)*bdp
    out.kc1=HArho(cHA*abc(3),[cHA abc])*bdp;
    out.bdp=bdp;

    % Find efficiency
    [~,eff]=HAc(cHA,v,'param',abc);
    out.eff=eff;

    if unonempty==true
        out.rho=HArho(u,c1);
        out.psi=HApsi(u,c1);
        out.psider=HApsider(u,c1);
        out.psix=HApsix(u,c1);
        out.wei=HAwei(u,c1);
    end


elseif strcmp(rhofunc,'HYP') || strcmp(rhofunc,'hyperbolic')
    % HYPERBOLIC
    out.class='HYP';

    param=[];
    if isempty(rhofuncparam)
        k=4.5;
    else
        k=rhofuncparam(1);
        if length(rhofuncparam) == 4
            param=rhofuncparam(2:4);
        end
    end

    if ~isempty(c)
        if isempty(param)
            [bdp, eff, A, B, d]=HYPc(c,v,'k',k);
        else
            [bdp, eff, A, B, d]=HYPc(c,v,'k',k, 'param', param);
        end
    elseif ~isempty(eff)
        [c, A, B, d]=HYPeff(eff,v,k);
        bdp=HYPc(c,v,'k',k, 'param', [A B d]);
    else
        [c, A, B, d]=HYPbdp(bdp,v,k);
        [~,eff]=HYPc(c,v,'k',k, 'param', [A B d]);
    end


    rhoHYPsup=HYPrho(200000,[c,k,A,B,d]);
    c1=[c,k,A,B,d];
    out.c1=c1;
    % kc = E(rho) = sup(rho)*bdp
    out.kc1=rhoHYPsup*bdp;
    out.bdp=bdp;
    out.eff=eff;

    if unonempty==true
        out.rho=HYPrho(u,c1);
        out.psi=HYPpsi(u,c1);
        out.psider=HYPpsider(u,c1);
        out.psix=HYPpsix(u,c1);
        out.wei=HYPwei(u,c1);
    end


elseif strcmp(rhofunc,'PD') || strcmp(rhofunc,'mdpd')
    % POWER DIVERGENCE

    if ~isempty(c)
        bdp=PDc(c);
    elseif ~isempty(eff)
        ceff=PDeff(eff);
        bdp=PDc(ceff);
    else
    end


    out.class='PD';
    c=PDbdp(bdp);
    out.c1=c;
    out.kc1=bdp;
    out.bdp=bdp;

    % Find efficiency
    [~,eff]=PDc(c);
    out.eff=eff;

    if unonempty==true
        out.rho=PDrho(u,c);
        out.psi=PDpsi(u,c);
        out.psider=PDpsider(u,c);
        out.psix=PDpsix(u,c);
        out.wei=PDwei(u,c);
    end

elseif strcmp(rhofunc,'HU') || strcmp(rhofunc,'huber')
    % HUBER

    if ~isempty(c)
        bdp=HUc(c,1);
    elseif ~isempty(eff)
        c=HUeff(eff,1);
        bdp=HUc(c,1);
    else
        if bdp==0.5
            c=0;
        else
            warning('FSDA:RhoPsiWei:WrongInputOpt','Huber link can have a bdp strictly greater than 0 only when c=0');
            c=NaN;
            bdp=NaN;
        end
    end


    out.class='HU';
    out.c1=c;
    out.kc1=Inf;
    out.bdp=bdp;

    % Find efficiency
    [~,eff]=HUc(c);
    out.eff=eff;

    if unonempty==true
        out.rho=HUrho(u,c);
        out.psi=HUpsi(u,c);
        out.psider=HUpsider(u,c);
        out.psix=HUpsix(u,c);
        out.wei=HUwei(u,c);
    end
else
    error('rho function not supported by code generation')
end


if unonempty==false
    out.rho=[];
    out.psi=[];
    out.psider=[];
    out.psix=[];
    out.wei=[];
end

end

%FScategory:UTISTAT