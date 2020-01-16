function out = boxtest(Y,group,varargin)
%boxtest performs Box test of equality of covariance matrices
%
%<a href="matlab: docsearchFS('boxtest')">Link to the help function</a>
%
% Box test is used to determine whether two or more covariance matrices
% are equal. Note that Box test is sensitive to departures from
% normality. If the samples come from non-normal distributions, then Box
% test may simply be testing for non-normality.
%
% Required input arguments:
%
% Y :           Input data. Matrix. n x v data matrix. n observations and v
%               variables. Rows of Y represent observations, and columns represent
%               variables. Missing values (NaN's) and infinite values
%               (Inf's) are allowed, since observations (rows) with missing
%               or infinite values will automatically be excluded from the
%               computations.
%               Data Types - single | double
%
% group :       Grouping variable.
%               Numeric vector | logical vector | character array | cell
%               array of character vectors.
%               group is a grouping variable defined as a categorical
%               variable, vector, character array, or cell array of
%               character vectors. Two observations are in the same group
%               if they have the same value in the group array. The
%               observations in each group represent a sample from a
%               population.
%
%
% Optional input arguments:
%
%     Fapprox : Test based on F approximation. Boolean. If Fapprox is
%               true, the asymptotic F distribution of the value of Box
%               test is also computed. On the other hand is
%               Fapprox is false (default) just the chi2 approximation is
%               computed.
%               Example - 'Fapprox',true
%               Data Types - logical
% dispresults : Display results. Boolean. If dispresults is
%               true, the value of the test and the associated p-value
%               will be shown on the screen. On the other hand is
%               dispresults is false (default) nothing is shown on the
%               screen.
%               Example - 'dispresults',true
%               Data Types - logical
%
%  Output:
%
%    out:   structure which contains the following fields
%           out.LR   =  scalar which contains Box test (uncorrected)
%                       for equality of covariances.
%                       This is -2ln M (see 'More About section' for the
%                       definition of M).
%  out.LRchi2approx  =  scalar which contains Box test (corrected)
%                       for equality of covariances. This version is
%                       called $\chi^2$ approximation of the Box test.
%                       This is $-2(1-c_1)\ln M$ (see further details for the
%                       definition of $c_1$ and $M$). This value must
%                       be compared with a $\chi^2$ with $0.5v(v+1)(g-1)$
%                       degrees of freedom.
% out.LRchi2approx_pval =  scalar which contains the p-value of
%                       $\chi^2$ approximation of Box test of homogeneity of
%                       covariances.
%     out.LRFapprox  =  scalar which contains the $F$ approximation of Box
%                       test of homogeneity of covariances. This field is
%                       given just if input option Fapprox is true.
%  out.LRFapprox_pval =  scalar which contains the p-value of
%                       $F$ approximation of Box test of homogeneity of
%                       covariances. This field is given just if input
%                       option Fapprox is true.
%       out.Spl       = pooled variance covariance matrix.
%
%
% More About:
%
% We assume independent samples of size $n_1$, $n_2$, $\ldots$,
% $n_g$ from a $v$-multivariate normal distribution. The hypothesis of
% equality of covariances is:
% 
% \[
% H_0= \Sigma_1 = \Sigma_2 = \ldots = \Sigma_g
% \]
%
% To make the test we calculate:
% 
% \[
% M=\frac{ |S_1|^{n_1-1} |S_2|^{n_2-1} \ldots |S_g|^{n_g-1} }{|S_{pl}|^{\sum_{i=1}^g n_i-1} }
% \]
%
% where $S_i$ is the covariance matrix of group $i$ and $S_{pl}$ is the
% pooled sample covariance matrix. It is clear that we must have $n_i-1>v$
% otherwise $|S_i|=0$ for some $i$ and $M$ would be zero. The statistic $M$
% is a modification of the likelihood ratio test and varies between 0 and 1
% with values near 1 favouring $H_0$ and values near 0 leading to the
% rejection of $H_0$ (see Rencher (2002) p. 256 for further details). The
% quantity $-2 \ln M$ is approximately distributed as a $\chi^2$
% distribution and is given in $\mbox{out.LR}$.  The quantity $-2(1-c_1) \ln  M$
% (where $c_1$ is a small sample correction factor) is usually called
% correct Box test and is approximately distributed as a $\chi^2$ with $0.5
% (g-1) v(v+1)$ degrees of freedom. We reject $H_0$ if 
% $$-2(1-c_1) \ln M  >\chi^2_{1-\alpha}$$ where $\chi^2_{1-\alpha}$ is the 
% $1-\alpha$ quantile. The value of this test is contained in
% $\mbox{out.LRchi2approx}$  and the corresponding p-value of this test is
% given in $\mbox{out.LRchi2approx_pval}$.
% Box also derived and F approximation to the test. This test is computed
% just if input option Fapprox is true. The value of this last test and the
% corresponding p-values are given in $\mbox{out.LRFapprox}$ and
% $\mbox{LRFapprox_pval}$.
% WARNING: if the absolute value of the determinant of the
% covariance matrix of any group is less than 1e-40,
% a missing value for LR test is reported.
%
%
% See also manova1.m, tkmeans.m
%
% References:
%
% Mardia, K.V., J.T. Kent, and J.M. Bibby (1979), "Multivariate Analysis,"
% Academic Press, London, pp. 140. [MKB].
%
% Rencher A.C. (2002), "Methods of Multivariate Analysis", 2nd edition,
% Wiley, New York, pp. 280-284.
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('boxtest')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Box test for Iris data.
    % load iris data
    load fisheriris
    % Compute Box test of equality of covariance matrices
    out=boxtest(meas,species);
%}

%{
    %% Box test for Iris data displaying results.
    % load iris data
    load fisheriris
    % Compute Box test of equality of covariance matrices
    out=boxtest(meas,species,'dispresults',true);
%}

%{
    %% Box test for Iris with option Fapprox.
    % load iris data
    load fisheriris
    % Compute Box test of equality of covariance matrices
    out=boxtest(meas,species,'dispresults',true,'Fapprox',true)
%}

%% Beginning of code
if nargin<1
    error('FSDA:boxtest:missingInputs','Input data matrix is missing')
end
if nargin<2
    error('FSDA:boxtest:missingInputs','Grouping variable is missing')
end

% test version for releases older than 2013b in order to use option upper inside cdf
vertest=verLessThan('matlab','8.2.0');

Fapprox=false;
dispresults= false;

if nargin > 2
    options=struct('Fapprox',Fapprox,'dispresults',dispresults);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:boxtest:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    Fapprox=options.Fapprox;
    dispresults=options.dispresults;
end


% Convert group to cell array from character array, make it a column
if (ischar(group))
    group = cellstr(group);
end

if (size(group, 1) == 1)
    group = group';
end

% Make sure inputs have the correct size
n = size(Y,1);
if (size(group,1) ~= n)
    error(message('FSDA:boxtest:InputSizeMismatch'));
end

% Remove missing Y columns first in case this eliminates a group
nonan = (sum(isnan(Y), 2) == 0);
Y = Y(nonan,:);
group = group(nonan,:);

% Convert group to indices 1,...,g and separate names
[groupnum, gnames] = grp2idx(group);
ngroups = length(gnames);

% Remove NaN values again
nonan = ~isnan(groupnum);
if (~all(nonan))
    groupnum = groupnum(nonan);
    Y = Y(nonan,:);
end

[n,p] = size(Y);
realgroups = ismember(1:ngroups,groupnum);
% g = number of groups
g = sum(realgroups);

% W = within sum of squares matrix
%  pooled sample covariance matrix Spl =Spl/(n-g)
W = zeros(p,p);

%sz = vector which will contain the size of the groups
sz=zeros(ngroups,1);

% Initialize scalar sVZ which will contain the quantity
% \sum_{i=1}^g (n_i-1) \ln | \S_i |
sVZ=0;
for j=1:ngroups
    r = find(groupnum == j);
    nr = length(r);
    sz(j)=nr;
    if (nr > 1)
        z = Y(r,:);
        xm = mean(z);
        z = z - xm(ones(nr,1),:);
        zz=z'*z;
        W = W + zz;
        %  VZi= unbiased variance covariance matrix for group i
        VZi=zz/(nr-1);
        dVZi= abs(det(VZi));
    else 
        dVZi = 0;
    end
    
    if dVZi < 1e-40
        out=NaN;
        return
    else
        sVZ= sVZ+(nr-1)*log(dVZi);
    end
end

% Spl = pooled sample covariance matrix
Spl=W/(n-g);

% Total sum of squares matrix would be computed as
% xm = mean(Y);
% Y = Y - xm(ones(n,1),:);
% T = Y'*Y;
% The Between sum of squares matrix is the difference
% B = T - W;

% LR= Likelihood ratio test statistic
% LR is -2 ln M in Rencher notation
LR= (n-g)*log(det(Spl))-sVZ;
%  gam = correction factor due to Box (1949)
%  see Mardia et al. (MKB) p. 140
%  or Rencher (2002) p. 282 eq. 7.20
% Remark:  -2 ln M for Rencher corresponds to eq. 5.3.27 of MKB
c1=(2*p^2+3*p-1)/(6*(p+1)*(g-1))*( sum(1./(sz-1)) -1/(n-g));
gam=1-c1;
LRchi2approx=LR*gam;
% LRchi2approx_pval = p value of the test
if vertest
    LRchi2approx_pval= 1-chi2cdf(LRchi2approx,0.5*p*(p+1)*(g-1));
else
    LRchi2approx_pval= chi2cdf(LRchi2approx,0.5*p*(p+1)*(g-1),'upper');
end

% Also compute if requested F approximation of the test and associated
% p-value
if Fapprox == true
    c2=((p-1)*(p+2)/(6*(g-1)))* ( sum(1./((sz-1).^2)) -1/(sum(sz-1))^2);
    a1=0.5*(g-1)*p*(p+1);
    a2=(a1+2)/abs(c2-c1^2);
    b1=(1-c1-a1/a2)/a1;
    b2=(1-c1+2/a2)/a2;
    
    if c2>c1^2
        LRFapprox=b1*LR;
    else
        LRFapprox=(a2*b2*LR)/(a1*(1-b2*LR));
    end
    if vertest
        LRFapprox_pval=1-fcdf(LRFapprox,a1,a2);
    else
        LRFapprox_pval=fcdf(LRFapprox,a1,a2,'upper');
    end
end

if dispresults == true
    disp('*****************************************************************')
    disp('Test of homogeneity of covariances (without correction factor)')
    disp(LR);
    disp('LRchi2approx = Chi2 approximation of Box test of homogeneity of covariances')
    disp(LRchi2approx)
    disp('p value of the Chi2 approximation of Box test')
    disp(LRchi2approx_pval)
    %   " This value must be compared with a X^2 r.v. with";
    %    v(v+1)(g-1)/2   degrees of freedom  */
    %   0.5*v*(v+1)*(g-1);; "degrees of freedom";
    if Fapprox== true
        disp('LRFapprox = F approximation of Box test of homogeneity of covariances')
        disp(LRFapprox)
        disp('p value of the F approximation of Box test')
        disp(LRFapprox_pval)
    end
    disp('*****************************************************************')
    
end

out=struct;
out.LR=LR;
% Store chi2 version of LR test
out.LRchi2approx=LRchi2approx;
out.LRchi2approx_pval=LRchi2approx_pval;
% Store F version of LR test
if Fapprox== true
    out.LRFapprox=LRFapprox;
    out.LRFapprox_pval=LRFapprox_pval;
end

out.Spl=Spl;
end
%FScategory:MULT-Multivariate