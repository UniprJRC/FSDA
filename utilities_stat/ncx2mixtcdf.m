function [qfval,varargout]= ncx2mixtcdf(c,n,lb,nc,varargin)
%ncx2mixtcdf cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))
%
%<a href="matlab: docsearchFS('ncx2mixtcdf')">Link to the help function</a>
%
%     given random variable $Q$ defined as
%
%     \[  
%     Q = \lambda_1 \chi^2_1 + \lambda_2 \chi_2 + ... + \lambda_k \chi_k +\sigma X_0
%     \]
%
%     where $\chi^2_1, ..., \chi^2_k$ are $k$ non central chi squared random variables,
%     with non centrality parameters $\delta_1, ..., \delta_k$ and degrees of
%     freedom $df_1, ..., df_k$.
%     and $X_0$ is a standard normal random variable, the purpose of this
%     routine is to compute $F_Q(x | df, delta) = P(Q < x)$ , that is the
%     cdf of $Q$ evaluated at $x$.
%     --------------------------------
%
% Required input arguments:
%
% c   :         value for which cdf has to be computed. Scalar. Value at
%               which the cdf must be evaluated
% n   :         degrees of freedom. Vector. Vector of length k containing
%               the degrees of freedom of the
%               k non central chi2 distributions
% lb  :         Coefficients of linear combination. Vector.
%               Vector of length k containing the coefficients of the
%               linear combinations of the k non central chi2 distributions
% nc  :         Non centrality parameters. Vector. Vector of length k
%               containing the k non centrality parameters
%               of the k non central chi2 distributions
%
%
% Optional input arguments:
%
% sigma :       standard deviation of N(0,1). Scalar. Coefficient
%               associated with standard deviation of the
%               standard normal distribution which can be added to the linear
%               combination of non central chi2 distributions
%               The default value of sigma is 0
%               Example - 'sigma',1
%               Data Types - double
%   lim :       Number of intergration terms. Scalar. Scalar which defines maximum number of integration terms.
%               The default value of lim is 10000
%               Example - 'lim',100000
%               Data Types - double
%   tol :       Tolerance. Scalar.
%               Scalar which controls the tolerance. The default value of
%               tolerance is 1e-09
%               Example - 'tol',1e-10
%               Data Types - double
%
% Remark:       The user should only give the input arguments that have to
%               change their default value.
%               The name of the input arguments needs to be followed by
%               their value. The order of the input arguments is of no
%               importance.
%
% Output:
%
%   qfval:      Value of cdf. Scalar. $qfval=F_Q(x | df, delta, sigma)$ is the value of the cdf of the mixture
%               evaluated at x
%
%  Optional Output:
%
%   tracert     : vector of length 7 containing
%                   tracert(1) = absolute sum
%                   tracert(2) = total number of integration terms
%                   tracert(3) = number of integrations
%                   tracert(4) = integration interval in final integration
%                   tracert(5) = truncation point in initial integration
%                   tracert(6) = standard deviation of initial convergence factor
%                   tracert(7) = number of iterations needed to locate integration parameters
%    ifault     :   scalar which informs about output of the procedure
%                   ifault=0 everything went OK
%                   ifault=1 required accuracy not achieved
%                   ifault=2 round-off error possibly significant
%                   ifault=3 invalid parameters (df or non centr parameters
%                               smaller than 0 or lmin=lmax=sigma=0)
%                   ifault=4 unable to locate integration parameters
%
%
%
% See also FSMenvmmd.m, FSM.m
%
% References:
%
%   Davies (1973), Numerical inversion of a characteristic function, vol.
%   60, Biometrika, pp. 415-417
%
%   Davies (1980), The distribution of a linear combination of Chi^2
%   Random variables, Applied Statistics vol. pp. 323-333
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('ncx2mixtcdf')">Link to the help function</a>
% Last modified 06-Feb-2015


% Examples:

%


%{
    % cdf of two chi squared rv.
    % Evaluate the cdf of the sum of two chi2 r.v. with degrees of freedom
    % 5 and 4 and coefficients of the linear combination 2 and 3 and non
    % centrality parameters 1 and 6
    x=35;
    df=[5;4];
    lb=[2;3];
    nc=[1;6];
    [out]=ncx2mixtcdf(x,df,lb,nc);
%}

%{
    % cdf of the sum of two non central chi2.
    % Evaluate the cdf of the sum of two chi2 r.v. with degrees of freedom
    % 5 and 4 and coefficients of the linear combination 2 and 3 and non
    % centrality parameters 1 and 6. Evaluate the cdf in a series of values
    % of x and plot the output
    df=[5;4];
    lb=[2;3];
    nc=[1;6];
    xx=0:1:100;
    cdfnc=zeros(length(xx),1);
    ij=1;
    for x=xx
        [out]=ncx2mixtcdf(x,df,lb,nc);
        cdfnc(ij)=out;
        ij=ij+1;
    end
    plot(xx',cdfnc)
    xlabel('x')
    ylabel('cdf of the mixture of non central X2')
%}

%{
    % Test tolerance.
    % Example which tests the results using different tolerances and
    % a different number of integration terms
    df=[1;1];
    lb=[-0.965785811006555;-0.681122597105154];
    nc=[0.2;0.3];
    x=-2.386488889335108;
    [out]=ncx2mixtcdf(x,df,lb,nc);
    disp('Value of cdf using default number of integration terms and default tolerance')
    disp(out)
    disp('-------------------------')
    tol=1e-06;
    [out]=ncx2mixtcdf(x,df,lb,nc,'tol',tol);
    disp(['Value of cdf using tol =' num2str(tol) ' and default integration terms'])
    disp(out)
    disp('-------------------------')
    lim=1000000;
    [out]=ncx2mixtcdf(x,df,lb,nc,'lim',lim);
    disp(['Value of cdf using numb. integration terms =' num2str(lim) ' and default tolerance'])
    disp(out)
    disp('-------------------------')
    lim=100000000;
    tol=1e-13;
    disp(['Value of cdf using numb. integration terms =' num2str(lim) ' and tolerance=' num2str(tol)])
    disp('In this last case it takes some seconds')
    [out]=ncx2mixtcdf(x,df,lb,nc, 'lim',lim, 'tol',tol);
    disp(out)
%}



%% Beginning of code

% Initialize tracert (vector which forms the additional optional output
tracert=zeros(7,1);

% Initialization of sigma (standard deviation of additional N(0,1))
sigmadef=0;

% Initialization of maximum number of integration terms
limdef=1e06;

% Initialization of tolerance
toldef=1e-8;

% store default values in the structure options
options=struct('sigma',sigmadef,'lim',limdef,'tol',toldef);

% check user options and update structure options
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:ncx2mixtcdf:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 4
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end

sigma=options.sigma;
xlim=options.lim;
acc=options.tol;
acc1=acc;

global intl;
global ersm;
intl=0;
ersm=0;

ifault=0;

%log28=log(2)/8;
log28=.0866;

%disp('----')
seqa=[3.696    3.360    3.080    2.800    2.640    2.400    2.200 ...
    2.000    1.848     1.680    1.540    1.400    1.320    1.200 1.1 1];

% nrep= replicates n, u times (in the columns)
nrep=repmat(n,1,length(seqa));

% Check input parameter values
if min(n) < 0  ||  min(nc) < 0
    return
end

lmax = max([lb;0]);
lmin= min([0;lb]);

if lmin == 0.0 && lmax == 0.0 && sigma == 0.0
    return
end

r=length(lb);

% th = vector which contains indexes of ordered elements of lb
[~,th]=sort(lb,'descend');

sigsq = sigma^2;

sd = sigsq + sum((lb.^2).*(2*n+4*nc));
mea= sum(lb.*(n+nc));

sd = sqrt(sd);

if lmax < -lmin
    almx= -lmin;
else
    almx=lmax;
end


% starting values for findu, ctff
utx = 16.0 / sd;
up = 4.5 / sd;
un = -up;

% Choose Delta such that
% max {Prob(X<x-2 \pi/Delta; X>x+2 \pi/delta}<0.5*acc1
% where acc1 is the integration error


if sd == 0.0
    if c>0
        qfval =1;
    else
        qfval =0;
    end
    return
end

% call subroutine findu in order to find truncation point
% Routine findu computes truncation point
% First time routine findu is called with no convergence factor
% truncation point with no convergence factor
utx=findu(utx, .5 * acc1);


% does convergence factor help */
if c ~= 0.0  && (almx > 0.07 * sd)
    
    [cfec,fail]=cfe(c);
    tausq = .25 * acc1 /cfec;
    
    if fail
        % fail = false ;
        % disp('Warning: non convergence in initial routine cfe')
    elseif truncation(utx, tausq) < .2 * acc1
        
        sigsq = sigsq + tausq;
        utx=findu(utx, .25 * acc);
        tracert(6) = sqrt(tausq);
    end
    
end

tracert(5) = utx; % truncation point of initial integration
acc1 = 0.5 * acc1;

% find RANGE of distribution, quit if outside this */

l1=1;

while l1
    
    [d1tmp,up] = ctff(acc1, up);
    d1=d1tmp- c;
    
    if (d1 < 0.0)
        qfval = 1.0;
        return
    end
    
    [d2tmp,un]=ctff(acc1, un);
    d2 = c - d2tmp;
    
    if (d2 < 0.0)
        qfval = 0.0;
        return
    end
    % find integration interval
    if d1>d2
        intv=2*pi/d1;
    else
        intv=2*pi/d2;
    end
    % calculate number of terms required for main and
    %   auxillary integrations
    xnt = utx / intv;
    % xntn = effective number of terms in the integration
    xntm = 3.0 / sqrt(acc1);
    
    if xnt > xntm * 1.5
        % parameters for auxillary integration
        
        if (xntm > xlim)
            % If number of terms of auxiliary integration is greater
            % than the maximum number of terms supplied by the user
            % than required accuracy cannot be achieved
            ifault = 1;
            qfval=NaN;
            varargout{2}=ifault;
            return
        end
        
        ntm = floor(xntm+0.5);
        % intv1 upper truncation point of integration divided by
        % number of terms
        intv1 = utx / ntm;
        x = 2.0 * pi / intv1;
        if (x <= abs(c))
            break
            % main integration
        end
        
        [cfe1,~]=cfe(c - x);
        [cfe2,fail2]=cfe(c + x);
        
        tausq = .33 * acc1 / (1.1 * (cfe1+cfe2));
        if fail2
            break
        end
        
        acc1 = .67 * acc1;
        % auxillary integration
        % ntm = number of terms = K
        % intv1 = Delta
        
        integrate(ntm, intv1, tausq, false);
        
        xlim = xlim - xntm;
        sigsq = sigsq + tausq;
        tracert(3) = tracert(3) + 1;
        tracert(2) = tracert(2) + ntm + 1;
        % find truncation point with new convergence factor */
        [utx]=findu(utx, .25 * acc1);
        acc1 = 0.75 * acc1;
        % goto l1;
        
    else
        l1=0;
    end
    
end

% main integration
% intv = integration interval in final integration
tracert(4) = intv;
if (xnt > xlim)
    disp('The number of supplied integration terms is too small for the required precision')
    disp('Please increase optional input parameter lim or relax the tolerance tol')
    qfval=nan;
    ifault = 1;
    varargout{2}=ifault;
    return
end

nt = floor(xnt+0.5);

integrate(nt, intv, 0, true);
tracert(3) = tracert(3) + 1;
tracert(2) = tracert(2) + nt + 1;
qfval = 0.5 - intl;
tracert(1) = ersm;

% 	 test whether round-off error could be significant
% 	   allow for radix 8 or 16 machines */
up=ersm; x = up + acc / 10.0;

if abs(x-up)<eps
    ifault=2;
end

tracert(7) = 0;


varargout{1}=tracert;

varargout{2}=ifault;

    function ut=findu(utx, accx)
        %utx finds u (upper truncation point in the integration)
        %such that truncation(u) < accx and truncation(u / 1.2) > accx
        %findu calls subroutine truncation
        
        
        ut = utx;
        u = ut / 4.0;
        if truncation(u, 0.0) > accx
            % if truncation(u)>accx then increase u up to when truncation(u)
            % becomes <= accx
            u=ut;
            while truncation(u, 0.0) > accx;
                
                ut = ut * 4.0;
                u=ut;
            end
            
        else
            % if truncation(u)<accx then decrease u up to when truncation(u)
            % still remains <accx
            ut=u;
            u=u/4;
            while (truncation(u,0)<=accx)
               
                ut=u;
                u=u/4;
            end
            
        end
        
        % The following loop has been replaced by a vectorized form of
        % truncation
        % divis has been replace by seqa
        % divis=[2.0,1.4,1.2,1.1];
        
        % This final loop is just to refine u
        %         utini=ut;
        %         for ii=1:4
        %             u = ut/divis(ii);
        %             if  truncation(u, 0.0)  <=  accx
        %                 ut=u;
        %                 %disp(ii)
        %             end
        %         end
        
        uchk=ut./seqa;
        tru=truncationv(uchk, 0.0);
        utchk=uchk(tru<accx);
        if ~isempty(utchk)
            ut=utchk(1);
%         else
%             ut=utini;
        end
    end


    function err=truncation(u, tausq)
        %truncation finds truncation error due to truncation at u
        % u can just be a scalar. The vectorized form of this function is
        % called truncationv. This function is called by findu
        
        sum2 = (sigsq + tausq) * u^2;
        
        prod1 = 2.0 * sum2;
        u = 2.0 * u;
        
        x = (u * lb).^2;
        sum1 = nc'* (x ./ (1.0 + x));
        xgt1=x>1;
        xxgt1=x(xgt1);
        nxgt1=n(xgt1);
        
        prod2 = sum(nxgt1.* log(xxgt1));
        
        log1x=log1(x, true );
        nlog1x=n.*log1x;
        prod3=sum(nlog1x(xgt1));
        prod1=prod1+ sum(nlog1x)-prod3;
        % Alternative statement to find prod1
        % xgt1n=~xgt1;
        % prod1=prod1+sum(nlog1x(xgt1n));
        
        % prod3 = nxgt1 * log1(xxgt1, true );
        s = sum(nxgt1);
        sum1 = 0.5 * sum1;
        prod2 = prod1 + prod2;  prod3 = prod1 + prod3;
        x = exp(-sum1 - 0.25 * prod2) / pi;
        y = exp(-sum1 - 0.25 * prod3) / pi;
        if s==0
            err1=1;
        else
            err1= x*2/s;
        end
        
        if prod3 >1
            err2 =2.5*y;
        else
            err2=1;
        end
        if (err2 < err1)
            err1 = err2;
        end
        x = 0.5 * sum2;
        
        if  ( x  <=  y )
            err2=1.0;
        else
            err2= y / x;
        end
        
        if err1<err2
            err=err1;
        else
            err=err2;
        end
    end

% vectorize function truncation
    function err=truncationv(u, tausq)
        %truncation finds truncation error due to truncation at u
        % u can be a scalar or a row vector
        % Accordingly err will be a scalar or a row vector
        % This function is called by findu
            
        sum2 = (sigsq + tausq) * u.^2;
        
        prod1 = 2.0 * sum2;
        u = 2.0 * u;
        
        x = (lb *u).^2;
        % sum1 = row vector
        sum1 = nc'* (x ./ (1.0 + x));
        xleq1=x(:)<=1;
        
        nlogx=nrep.*log(x);
        % Slower alternative
        % nlogx=bsxfun(@times,n,log(x));
        
        nlogx(xleq1)=0;
        
        prod2 = sum(nlogx,1);
        
        log1x=log1(x, true );
        
        nlog1x=nrep.*log1x;
        % Slower alternative
        % nlog1x=bsxfun(@times,n,log1x);
        
        
        sumall=sum(nlog1x,1);
        nlog1x(xleq1)=0;
        prod3=sum(nlog1x,1);
        prod1=prod1+ sumall-prod3;
        
        % prod3 = nxgt1 * log1(xxgt1, true );
        
        nrep(xleq1)=0;
        s = sum(nrep,1);

        sum1 = 0.5 * sum1;
        prod2 = prod1 + prod2;  prod3 = prod1 + prod3;
        x = exp(-sum1 - 0.25 * prod2) / pi;
        y = exp(-sum1 - 0.25 * prod3) / pi;
        
        err1= x*2./s;
        err1(s==0)=1;
        
        err2 =2.5*y;
        err2(prod3<=1)=1;
        
        err1f=err1;
        err21=err2 < err1;
        err1f(err21)=err2(err21);
        
        x = 0.5 * sum2;
        
        
        err2f= y./ x;
        err2f(x<=y)=1;
        
        err=err2f;
        err12f=err1f<err2f;
        err(err12f)=err1f(err12f);
    end

    function s=log1(x,first)
        % x can be a scalar a vector or a matrix. Output s accordingly will
        % be a scalar a vector or a matrix 
        % In other words log1 function operates elementwise on arrays
        %  for the elements of x which are greater than 0.1
        % if (first) s= log(1 + x) ; else
        % s= log(1 + x) - x 
        
        [nx,cx]=size(x);
        % If x is a row vector than x is preliminarly transformed into a
        % column vector. Without this if the instruction below
        % s1=sabsxn+term/k; which appear before the while loop generates an
        % error because sabsxn is a column vector and on the other hand
        % term is a row vector.
        
        if nx==1 && cx >1
            x=x';
        end
        
        absx=abs(x(:))>0.1;
        
        xabsx=x(absx);
        s=zeros(nx*cx,1);
        
        absxn=~absx;
        xabsxn=x(absxn);
        y = xabsxn./(2 + xabsxn);
        
        term=2*exp(3*log(y));
        % Option below to compute term seems slower
        % term = 2 * y.^3;
        
        if first
            s(absx)=log(1.0 + xabsx);
            s(absxn) = 2*y;
        else
            s(absx)=log(1.0 + xabsx) - xabsx;
            s(absxn)= -xabsxn.*y;
        end
        
        if ~isempty(y) 
            k = 3;
            y = y.^2;
            
            sabsxn=s(absxn);
            s1=sabsxn+term/k;
            
            while   ~isequal(s1,sabsxn)  % max(abs(s1-sabsxn))~=0 %
                k = k + 2.0;
                term = term.* y;
                sabsxn = s1;
                s1=sabsxn+term/k;
            end
            s(absxn)=s1;
        end
        
        if cx>1
            s=reshape(s,nx,cx);
        end
        
    end

    function [coefftau,fail]=cfe(x)
        %  coef of tausq in error when convergence factor of
        %   exp1(-0.5*tausq*u^2) is used when df is evaluated at x
        % This function computes eq (10) (integration error)
        % fail =boolean (true if integration error is greater then 100)
        
        axl = abs(x);
        
        if x >0
            sxl=1;
        else
            sxl=-1;
        end
        sum1 = 0.0; % row 224
        
        for  jj =r:-1:1
            t = th(jj);
            if lb(t) * sxl > 0.0 
                lj = abs(lb(t));
                axl1 = axl - lj * (n(t) + nc(t));
                axl2 = lj / log28;
                
                if axl1 > axl2 
                    axl = axl1;
                else
                    if axl > axl2
                        axl = axl2;
                    end
                    
                    sum1 = (axl - axl1) / lj;
                    
                    thj1=th(1:jj-1);
                    sum1=sum1+sum(n(thj1)+nc(thj1));
                    
                    break
                end
            end
            
        end
        if (sum1 > 100.0)
            fail =true;
            coefftau=  1.0;
        else
            coefftau=2^(sum1 / 4.0) / (pi * (axl^2));
            fail =false;
        end
    end


    function [c2,upn] =ctff(accx,upn)
        %ctff find cut off (say c2)  so that p(qf > c2) < accx  if (upn > 0,
        %  p(qf < c2) < accx otherwise */
        
        u2 = upn;   u1 = 0.0;  c1 = mea;
        if u2>0
            rb=2*lmax;
        else
            rb=2*lmin;
        end
        
        u=u2/(1+u2*rb);
        [out,c2]=errbd(u);
        while(out>accx)
            u1=u2;
            c1=c2;
            u2=2*u2;
            
            u = u2 / (1.0 + u2 * rb);
            [out,c2]=errbd(u);
        end
        
        
        u=(c1 - mea) / (c2 - mea);
        while u < 0.9
            u=(u1+u2)/2;
            
            [out,xconst]=errbd(u/(1+u*rb));
            if out>accx
                u1=u;
                c1=xconst;
            else
                u2=u;
                c2=xconst;
            end
            u=(c1 - mea) / (c2 - mea);
        end
        upn=u2;
    end

    function [out,cx]=errbd(u)
        %  find bound on tail probability using mgf, cutoff
        %  point returned to *cx */
        xconst = u * sigsq;
        sum1 = u * xconst;
        u = 2.0 * u;
        
        x = u * lb; y = 1.0 - x;
        cx = xconst+lb' * ((nc./y + n) ./ y);
        x2=x.^2;
        sum1 =  sum1+nc'* (x2./(y.^2)) + n' * (x2./ y + log1(-x, false));
        
        
        % Alternative way to find sum1 (seems slower)
        % sum1 =  u * xconst+ sum(nc.* (x2./(y.^2)) + n.* (x2./ y + log1(-x, false)));
        % sum1 =  sum1+nc'* ((x./y).^2) + n' * (x.^2./ y + log1(-x, false));

        % Also this alternative to find sum1 seems slower
        % x2y=x2./y.^2;
        % sum1chk =  sum1chk+nc'* (x2y) + n' * (x2y.*y + log1(-x, false));

        out=exp(-0.5 * sum1);
    end


    function integrate(nterm, interv, tausq, mainx)
        %integrate computes auxiliary and final integration
        %integrate carries out integration with nterm terms at stepsize
        %interv
        %
        %  Required input arguments:
        %
        %    nterm = scalar. Number of terms (K)
        %    interv = scalar: stepsize (Delta)
        %    tausq = scalar which represents the standard deviation of the
        %    additional normal variable \tau Z which is added to the sum of
        %    mixture of non central chi^2. This acts as a sort of
        %    convergence factor and it can help to reduce the truncation
        %    point U
        %    mainx is a boolean. if mainx is true the integrand is
        %    multiplied by 1.0-exp(-0.5*tausq*u^2)
        %
        %  Output:
        %   intl = output of the integration
        %        intl = \sum_1^K Im[ \phi \{ u e^-i u*c \}]/(\pi \Delta u)
        %    u=(k+0.5)\Delta
        %
        %  Espression \sum_1^K Im[ \phi \{ u e^-i u*c \}] is equal to the
        %  exp { } expression in curly brackets
        
        
        inpi = interv / pi;
        k = nterm:-1:0;
        u = (k + 0.5) * interv;
        sum1 = - 2.0 * u * c;
        sum2 = abs(sum1);
        sum3 = - 0.5 * sigsq * u.^2;
        
        x = 2.0 * lb * u;  y = x.^2;
        sum3 = sum3 - 0.25 * n'* log1(y, true);
        
        
        y= bsxfun(@times,nc,x)./ (1.0 + y);
        z =  bsxfun(@times,n,atan(x))+y;
        
        %             y = nc.* x ./ (1.0 + y);
        %             z = n.* atan(x) + y;
        sum1 = sum1 + sum(z,1);   sum2 = sum2 + sum(abs(z),1);
        sum3 = sum3 - 0.5 *sum(x.*y,1);
        
        
        x = inpi * exp(sum3) ./ u;
        if ~mainx
            x = x.* (1.0 - exp(-0.5 * tausq * u.^2 ));
        end
        sum1 = sin(0.5 * sum1) .* x;
        sum2 = 0.5 * sum2 .* x;
        intl = intl + sum(sum1);
        ersm = ersm + sum(sum2);
        
    end
end