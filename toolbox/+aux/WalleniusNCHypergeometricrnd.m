function r = WalleniusNCHypergeometricrnd(n,m,N,odds, accuracy)
%WalleniusNCHypergeometricrnd returns a random number from non-central hypergeometric probability density function
% n = number of balls taken
% m = number of red balls
% N = total number of balls in the urn
% omega = odds
if nargin<4
    odds=1;
    accuracy=1e-08;
elseif  nargin<5
    accuracy=1e-08;
else

end
cutoff=1e-10;


if (n >= N || m >= N || n <= 0 || m <= 0 || odds <= 0.)
    % // trivial cases
    if (n == 0 || m == 0)
        r=0;
    elseif (m == N)
        r=n;
    elseif n == N
        r=m;
    else
    end

    return
end

if (odds == 0.)
    if (n > N-m)
        error('FSDA:WalleniusNCHypergeometricrnd:WrgInpt',"Not enough items with nonzero weight in function WalleniusNCHyp");
    end
end


if (odds == 1.)
    % // use hypergeometric function if odds == 1
    % M. K. N
    r = hygernd(N, m, n);
    return
end


if (n < 30)
    r=WalleniusNCHypUrn(n, m, N, odds);

elseif (n*N < 10000)
    r=WalleniusNCHypTable(n, m, N, odds, cutoff, accuracy);

else
    r= WalleniusNCHypRatioOfUnifoms(n, m, N, odds, accuracy);
end
end


%// sampling from Wallenius noncentral hypergeometric distribution
%// by simulating urn model
function x=WalleniusNCHypUrn (n, m, N, odds)
%int32_t x;                           // sample
%int32_t m2;                          // items of color 2 in urn
%double mw1, mw2;                     // total weight of balls of color 1 or 2
x = 0;
m2 = N - m;
mw1 = m * odds;
mw2 = m2;
while n>0
    randFromU01=unifrnd(0,1);
    if (randFromU01 * (mw1 + mw2) < mw1)
        x=x+1;
        m=m-1;
        if (m == 0)
            break;
        else
            mw1 = m * odds;
        end
    else
        m2=m2-1;
        if (m2 == 0)
            x = x+n-1;
            break
        end
        mw2 = m2;
    end
    n=n-1;
end
end




function x=WalleniusNCHypTable(n, m, N, odds, cutoff, accuracy)
% // Sampling from Wallenius noncentral hypergeometric distribution
% // using chop-down search from a table created by recursive calculation.
% // This method is fast when n is low or when called repeatedly with
% // the same parameters.


% if (n != wnc_n_last || m != wnc_m_last || N != wnc_N_last || odds != wnc_o_last) {
%    // set-up: This is done only when parameters have changed
%    wnc_n_last = n;  wnc_m_last = m;  wnc_N_last = N;  wnc_o_last = odds;

% CWalleniusNCHypergeometric wnch(n,m,N,odds);   // make object for calculation
% xfirst,xlast
[success, wall_ytable, wall_x1, x2] = MakeTable(n,m,N,odds, cutoff); %  wall_ytable, WALL_TABLELENGTH); % // make table of probability values

if success == true
    wall_tablen = x2 - wall_x1 + 1;         % // table long enough. remember length
else
    wall_tablen = 0;                      %   // remember failure
end


if (wall_tablen == 0)
    % // table not long enough. Use another method
    x=WalleniusNCHypRatioOfUnifoms(n,m,N,odds, accuracy);
end

while 1>0                             %    // repeat in the rare case of failure
    u = unifrnd(0,1);                              %// uniform variate to convert
    for xx=1:wall_tablen           % // chop-down search
        u= u- wall_ytable(xx);
        if (u < 0)
            x= xx + wall_x1-1;
            return       % // value found
        end
    end
end
end

function [i1returned, tabl, xfirst, xlast]=MakeTable(n,m,N,omega,cutoff)
% // Makes a table of Wallenius noncentral hypergeometric probabilities
%   // table must point to an array of length MaxLength.
%   // The function returns 1 if table is long enough. Otherwise it fills
%   // the table with as many correct values as possible and returns 0.
%   // The tails are cut off where the values are < cutoff, so that
%   // *xfirst may be > xmin and *xlast may be < xmax.
%   // The value of cutoff will be 0.01 * accuracy if not specified.
%   // The first and last x value represented in the table are returned in
%   // *xfirst and *xlast. The resulting probability values are returned in
%   // the first (*xfirst - *xlast + 1) positions of table. Any unused part
%   // of table may be overwritten with garbage.
%   //
%   // The function will return the following information when MaxLength = 0:
%   // The return value is the desired length of table.
%   // *xfirst is 1 if it will be more efficient to call MakeTable than to call
%   // probability repeatedly, even if only some of the table values are needed.
%   // *xfirst is 0 if it is more efficient to call probability repeatedly.

% double * p1, * p2;                  // offset into p
% double mxo;                         // (m-x)*omega
% double Nmnx;                        // N-m-nu+x
% double y, y1;                       // probability. Save old p[x] before it is overwritten
% double d1, d2, dcom;                // divisors in probability formula
% double area;                        // estimate of area needed for recursion method
% int32_t xi, nu;                     // xi, nu = recursion values of x, n
% int32_t x1, x2;                     // lowest and highest x or xi
% int32_t i1, i2;                     // index into table
% int32_t UseTable;                   // 1 if table method used
% int32_t LengthNeeded;               // Necessary table length
%
% // special cases
% if (n == 0 || m == 0) {x1 = 0; goto DETERMINISTIC;}
% if (n == N)           {x1 = m; goto DETERMINISTIC;}
% if (m == N)           {x1 = n; goto DETERMINISTIC;}
% if (omega <= 0.) {
%    if (n > N-m) FatalError("Not enough items with nonzero weight in  CWalleniusNCHypergeometric::MakeTable");
%    x1 = 0;
%    DETERMINISTIC:
%    if (MaxLength == 0) {
%       if (xfirst) *xfirst = 1;
%       return 1;
%    }
%    *xfirst = *xlast = x1;
%    *table = 1.;
%    return 1;
% }

xmin = m + n - N;
% calculate xmin
if (xmin < 0)
    xmin = 0;
end
xmax = n;
%  calculate xmax
if (xmax > m)
    xmax = m;
end

accuracy=1e-08;
MaxLength=512;

if (cutoff <= 0. || cutoff > 0.1)
    cutoff = 0.01 * accuracy;
end

LengthNeeded = N - m;               % // m2
if (m < LengthNeeded)
    LengthNeeded = m;
end

if (n < LengthNeeded)
    LengthNeeded = n; % // LengthNeeded = min(m1,m2,n)
end


% if MaxLength <= 0
%     % // Return UseTable and LengthNeeded
%     if xfirst *xfirst == UseTable
%         i1 = LengthNeeded + 2;           % Necessary table length
%         if ~UseTable && i1 > 200
%             % // Calculate necessary table length from standard deviation
%             % double sd = sqrt(variance()); // calculate approximate standard deviation
%             % // estimate number of standard deviations to include from normal distribution
%             i2 = NumSD(accuracy) * sd + 0.5;
%             if (i1 > i2)
%                 i1 = i2;
%             end
%             i1returned=i1;
%             return
%         end
%     end
% end

tabl=zeros(LengthNeeded,1);
p1=tabl;
p2=p1;

% area = n*LengthNeeded;     % // Estimate calculation time for table method
% UseTable = area < 5000 || (area < 10000. && N > 1000. * n);

% if (UseTable && MaxLength > LengthNeeded)
    % // use recursion table method
    p1(1) = 1;
    p2(1)=p1(1);
    tabl(1)=0;
    tabl(2)=1;

    %  // make space for p1[-1]
    % p1[-1] = 0.;  p1[0] = 1.;        // initialize for recursion
    x1 = 0;
    x2 = 0;
    for nu = 1:n
        if (n - nu < xmin - x1 || p1(x1+1) < cutoff)
            x1=x1+1;                     % // increase lower limit when breakpoint passed or probability negligible
            p2=[0;p2];                     %#ok<AGROW> % // compensate buffer offset in order to reduce storage space
            remove1=1;
        else
            remove1=0;
        end
        if (x2 < xmax && p1(x2+1) >= cutoff)
            x2=x2+1;
            y1 = 0.;           %  // increase upper limit until x has been reached
        else
            y1 = p1(x2+1);
        end
        if (p2(1) - tabl(1) + x2 >= MaxLength || x1 > x2)
                error('FSDA:WalleniusNCHypergeometricrnd','table length exceeded') % goto ONE_BY_ONE;       %     // Error: table length exceeded. Use other method
        end

        mxo = (m-x2)*omega;
        Nmnx = N-m-nu+x2+1;
        for xi = x2:-1:x1
            % // backwards loop
            d2 = mxo + Nmnx;
            mxo = mxo+ omega;
            Nmnx= Nmnx-1;
            d1 = mxo + Nmnx;
            dcom = 1. / (d1 * d2);     % // save a division by making common divisor
            if xi>0
                y  = p1(xi)*mxo*d2*dcom + y1*(Nmnx+1)*d1*dcom;
                y1 = p1(xi);             % // (warning: pointer alias, can't swap instruction order)
            else
                y=y1*(Nmnx+1)*d1*dcom;
                y1=0;
            end
            p2(xi+1) = y;
            % try
            % tabl(2:end)= p2(1:end-1);
            % catch
            %     ddd=1;
            % end
            % tabl(xi+3)
            p1(xi+1-remove1)=y;
            tabl=p1;
            %disp(['nu='  num2str(nu)])
            %disp(tabl)
        end
        p1 = p2;
    end

    tabl=p2;
    % // return results
    i1 =  x2 - x1 + 1;
    i2 = i1; %// desired table length
    if (i2 > MaxLength)
        i2 = MaxLength; % // limit table length
    end

    xfirst = x1;
    xlast = x1 + i2 - 1;

    % if (i2 > 0)
    %       tabl=tabl(i2+1:end); %// copy to start of table
    %  end

    ini=find(tabl~=0,1,'first');
    tabl=tabl(ini:end,1);

    i1returned=i1 == i2;                  %  // true if table size not reduced
% else
%     % // Recursion method would take too much time
%     %// Calculate values one by one
%     % ONE_BY_ONE:
% 
%     % // Start to fill table from the end and down. start with x = floor(mean)
%     % x2 = (int32_t)mean();
%     x1 = x2 + 1;
%     i1 = MaxLength;
%     while (x1 > xmin)               % // loop for left tail
%         x1=x1-1;  i1=i1-1;
%         y = WalleniusNCHypergeometricpdf(x1,n,m,N,odds,accuracy);
%         tabl(i1) = y;
%         if (y < cutoff)
%             break
%         end
% 
%         if (i1 == 0)
%             break
%         end
% 
%         xfirst = xfirst*x1;
%         i2 = x2 - x1 + 1;
%         if (i1 > 0 && i2 > 0)
%             % // move numbers down to beginning of table
%             memmove(table, table+i1, i2*sizeof(table(0)));
%         end
%         % // Fill rest of table from mean and up
%         i2=i2-1;
%         while (x2 < xmax)            %    // loop for right tail
%             if (i2 == MaxLength-1)
%                 xlast = xlast*x2;
%                 i1returned=0;
%                 return     % // table full
%             end
%             x2=x2+1;  i2=i2+1;
%             y = WalleniusNCHypergeometricpdf(x2,n,m,N,odds,accuracy);
%             tabl(i2) = y;
%             if (y < cutoff)
%                 break
%             end
%             xlast = xlast*x2;
%             i1returned=1;
%             return
%         end
%     end
% end
end

function xi=WalleniusNCHypRatioOfUnifoms(n,  m,  N, odds, accuracy)
%    // sampling from Wallenius noncentral hypergeometric distribution
%    // using ratio-of-uniforms rejection method.

%
%    // Make object for calculating mean and probability.
%    CWalleniusNCHypergeometric wnch(n, m, N, odds, accuracy);
%
xmin = m+n-N;
if (xmin < 0)
    xmin = 0;  % // calculate limits
end

xmax = n;
if (xmax > m)
    xmax = m;
end

mea=CWalleniusNCHypergeometricmean(n,m,N,odds,xmin,xmax);
%
%       // find approximate mean
%       mean = wnch.mean();
%
% // find approximate variance from Fisher's noncentral hypergeometric approximation
r1 = mea * (m-mea); r2 = (n-mea)*(mea+N-n-m);
variance = N*r1*r2/((N-1)*(m*r2+(N-m)*r1));


UseChopDown = variance < 4;       % // use chop-down method if variance is low

if UseChopDown == false
    % // find mode (same code in CWalleniusNCHypergeometric::mode)
    wnc_mode=floor(mea);
    f2 = 0.;
    if (odds < 1)
        if (wnc_mode < xmax)
            wnc_mode= wnc_mode+1;
        end

        x2 = xmin;
        if (odds > 0.294 && N <= 10000000)
            x2 = wnc_mode - 1;
        end

        % // search for mode can be limited
        for xi = wnc_mode:-1:x2
            f = aux.WalleniusNCHypergeometricpdf(xi,n,m,N,odds,accuracy);
            if (f <= f2)
                break
            end
            wnc_mode = xi;
            f2 = f;
        end
    else
        if (wnc_mode < xmin)
            wnc_mode=wnc_mode+1;
        end

        x2 = xmax;
        if (odds < 3.4 && N <= 10000000)
            x2 = wnc_mode + 1;                %     // search for mode can be limited
        end

        for xi = wnc_mode:1:x2
            f = aux.WalleniusNCHypergeometricpdf(xi,n,m,N,odds,accuracy);
            if (f <= f2)
                break
            end

            wnc_mode = xi; f2 = f;
        end
    end
    wnc_k = f2; %  value at mode

    % // find approximate variance from normal distribution approximation
    rsqrt2pi=1/sqrt(2*pi);
    variance = rsqrt2pi / wnc_k;
    variance = variance^2;

    % find center and width of hat function
    wnc_a = mea + 0.5;
    s123 = 0.40 + 0.8579*sqrt(variance+0.5) + 0.4*abs(mea-wnc_mode);
    s4 = 0.;
    r1 = xmax - mea - s123;  r2 = mea - s123 - xmin;
    if (r1 > r2)
        r1 = r2;
    end

    if ((odds>5. || odds<0.2) && r1>=-0.5 && r1<=8)
        % // s4 correction needed
        if r1 < 1
            r1 = 1;
        end
        s4 = 0.029 * N^0.23 / (r1*r1);
    end

    wnc_h = 2. * (s123 + s4);

    %// find safety bounds
    wnc_bound1 = floor(mea - 4. * wnc_h);
    if wnc_bound1 < xmin
        wnc_bound1 = xmin;
    end

    wnc_bound2 = floor(mea + 4. * wnc_h);
    if (wnc_bound2 > xmax)
        wnc_bound2 = xmax;
    end

else  % UseChopDown == true
    accuracy=1e-9;
    % // for small variance, use chop down inversion
    xi= WalleniusNCHypInversion(n,m,N,odds, accuracy,mea);
    return
end

% // use ratio-of-uniforms rejection method
while 1>0                                    %  rejection loop
    u = unifrnd(0,1);
    if (u == 0.)
        continue                    % // avoid division by 0
    end

    uagain=unifrnd(0,1);
    x = wnc_a + wnc_h * (uagain-0.5)/u;
    if (x < 0. || x > 2E9)
        continue           % // reject, avoid overflow
    end

    xi = floor(x);                       %   // truncate
    if (xi < wnc_bound1 || xi > wnc_bound2)
        continue                               % // reject if outside safety bounds
    end
    % #if 0 // use rejection in x-domain
    % // use rejection in t-domain (this is faster)
    % double hx, s2, xma2;                       // compute h(x)
    s2 = wnc_h * 0.5;
    s2 = s2^2;
    xma2 = xi - (wnc_a-0.5);
    xma2 = xma2^2;

    if (s2 >= xma2)

        hx = 1;
    else
        hx= s2 / xma2;
    end

    % // rejection in t-domain implemented in CWalleniusNCHypergeometric::BernouilliH
    r=wnchBernouilliH(hx * wnc_k * 1.01, u * u * wnc_k  * 1.01, n,xi,N,m, odds);
    if r == true
        break                             %     // acceptance
    end
    %     // rejection
end
end


function mea=CWalleniusNCHypergeometricmean(n,m,N,omega,xmin,xmax)
if omega == 1
    % { // simple hypergeometric
    mea= m*n/N;
    return
end

if (omega == 0.)
    if  (n > N-m)
        error('FSDA:WalleniusNCHypergeometricrnd:WrgInp',"Not enough items with nonzero weight in CWalleniusNCHypergeometric::mean");
    end
end

if (xmin == xmax)
    mea=xmin;
    return
end

% // calculate Cornfield mean of Fisher noncentral hypergeometric distribution as first guess
a = (m+n)*omega + (N-m-n);
b = a*a - 4.*omega*(omega-1.)*m*n;
if b>0
    b = sqrt(b);
else
    b=0;
end

mea = (a-b)/(2.*(omega-1.));
if (mea < xmin)
    mea = xmin;
end

if (mea > xmax)
    mea = xmax;
end

m1r = 1./m;  m2r = 1./(N-m);
iter = 0;

mea1=Inf;
if (omega > 1.)
    while abs(mea1 - mea) > 2E-6

        % // Newton Raphson iteration
        mea1 = mea;
        e1 = 1.-(n-mea)*m2r;
        if e1 < 1E-14
            e2 = 0.;     % // avoid underflow
        else
            e2 = e1^(omega-1);
        end

        g = e2*e1 + (mea-m)*m1r;
        gd = e2*omega*m2r + m1r;
        mea = mea -g / gd;
        if (mea < xmin)
            mea = xmin;
        end

        if (mea > xmax)
            mea = xmax;
        end
        iter=iter+1;

        if iter > 40
            error('FSDA:WalleniusNCHypergeometricrnd:WrhIter',"Search for mean failed in function CWalleniusNCHypergeometric::mean");
        end
    end

else % { // omega < 1
    omegar = 1./omega;

    while abs(mea1 - mea) > 2E-6

        mea1 = mea;
        e1 = 1.-mea*m1r;
        if (e1 < 1E-14)
            e2 = 0.;   % // avoid underflow
        else
            e2 = e1^(omegar-1);
        end
        g = 1.-(n-mea)*m2r-e2*e1;
        gd = e2*omegar*m1r + m2r;
        mea = mea-g / gd;
        if (mea < xmin)
            mea = xmin;
        end

        if (mea > xmax)
            mea = xmax;
        end
        iter=iter+1;
        if (iter > 40)
            error('FSDA:WalleniusNCHypergeometricrnd:WrngIte',"Search for mean failed in function CWalleniusNCHypergeometric::mean");
        end
    end
end
end


function r=WalleniusNCHypInversion(n, m, N, odds, accuracy, mea)
% // sampling from Wallenius noncentral hypergeometric distribution
% // using down-up search starting at the mean using the chop-down technique.
% // This method is faster than the rejection method when the variance is low.

% Make objects for calculating mean and probability.
% It is more efficient to have two identical objects, one for down search
% and one for up search, because they are obtimized for consecutive x values.
% CWalleniusNCHypergeometric wnch1(n, m, N, odds, accuracy);
% CWalleniusNCHypergeometric wnch2(n, m, N, odds, accuracy);

accura = accuracy * 0.01;
if (accura > 1E-7)
    accura = 1E-7;            %  // absolute accuracy
end

wall_x1 = floor(mea);            % // start at floor and ceiling of mean
x2 = wall_x1 + 1;
xmin = m+n-N;
if (xmin<0)
    xmin = 0;           % // calculate limits
end

xmax = n;
if (xmax>m)
    xmax = m;
end

updown = 3;                                   %// start searching both up and down
% ins=[2:4:30 3:4:31];
ins=[2 3];
while 1>0                                    %// loop until accepted (normally executes only once)
    u = unifrnd(0,1);                           %    // uniform random number to be converted
    while (updown>0)                           % // search loop
        if mod(updown,2)==1 % updown is odd     %                   // search down
            if (wall_x1 < xmin)
                % if mod(updown,2)==1
                    updown=updown-1; % // stop searching down
                % end
            else
                f = aux.WalleniusNCHypergeometricpdf(wall_x1,n,m,N,odds,accuracy);
                u = u-f;                           % // subtract probability until 0
                if (u <= 0)
                    r=wall_x1;
                    return
                else
                    wall_x1=wall_x1-1;
                end
                if (f < accura)
                    % if mod(updown,2)==1
                        updown = updown-1;     %% // stop searching down
                    %end
                end
            end
        end

        % updown =2,3,6,7,10,11,... %  (updown & 2) 
        if sum(updown==ins)==1         %   // search up
            if x2 > xmax
                % if sum(updown==ins)==1
                    updown =updown-2;        %  updown &= ~2;  // stop searching up
                % end
            else
                f = aux.WalleniusNCHypergeometricpdf(x2,n,m,N,odds,accuracy);
                u = u-f;                           % // subtract probability until 0
                if u <= 0
                    r=x2;
                    return
                else
                    x2=x2+1;
                    if (f < accura)
                        %if sum(updown==ins)==1
                            updown=updown-2;    %  // stop searching down
                        % end
                    end
                end
            end
        end
    end
end
end


function testnum=wnchBernouilliH(h,rh, n,x_,N,m,odds)
% // This function generates a Bernouilli variate with probability proportional
%    // to the univariate Wallenius' noncentral hypergeometric distribution.
%    // The return value will be 1 with probability f(x_)/h and 0 with probability
%    // 1-f(x_)/h.
%    // This is equivalent to calling sto->Bernouilli(probability(x_)/h),
%    // but this method is faster. The method used here avoids calculating the
%    // Wallenius probability by sampling in the t-domain.
%    // rh is a uniform random number in the interval 0 <= rh < h. The function
%    // uses additional random numbers generated from sto.
%    // This function is intended for use in rejection methods for sampling from
%    // the Wallenius distribution. It is called from
%    // StochasticLib3::WalleniusNCHypRatioOfUnifoms in the file stoc3.cpp

x = x_;                         % // save x in class object
bico=lnbico(n,x,N,m);
% // calculate bico = log(Lambda)
[rd,E,r,w]=findpars(n,x,N,m,odds);
omegai=[odds,1];
xi=[x, n-x];

if (E > 0.)
    k = log(E);           %       // correction for majorizing function
    k = 1. + 0.0271 * (k * sqrt(k));
else
    k = 1.;
end
k = k*w;                    %     // w * k
rdm1 = rd - 1;

% // calculate phi(½)/rd
LN2=log(2);
phideri0 = -LN2 * rdm1;
for i=1:2
    romegi = r * omegai(i);
    if romegi > 40
        qi=0;
         qi1 = 1;           % // avoid underflow
    else
        [qi1,qi] = pow2_1(-romegi);
    end
    phideri0 = phideri0+ xi(i) * log1mx(qi, qi1);
end

rsqrt8=1/sqrt(8);
sqrt2pi=sqrt(2*pi);
erfk = erf(rsqrt8 / k);
f0 = rd * exp(phideri0 + bico);
G_integral = f0 * sqrt2pi * k * erfk;

ts=Inf;
if (G_integral <= h)          % // %G fits under h-hat
    while (abs(ts) >= 0.5)
        ts = normrnd(0,k,1,1);    %   // sample ts from normal distribution
       % ts= -0.015890922832026869
        %  // reject values outside interval, and avoid ts = 0
    end

        ts = ts+0.5;                 %   // ts = normal distributed in interval (0,1)

    fts=0;
    for j=1:2
        % // calculate (Phi(ts)+Phi(1-ts))/2
        logts = log(ts);
        rlogts = r * logts; % // (ts = 0 avoided above)

        % log1pow0=log((1-exp(rlogts*odds))^xi(1));
        % log1pow1=log((1-exp(rlogts))^xi(2));
        log1pow0=log1pow(rlogts*odds,xi(1));
        log1pow1=log1pow(rlogts,xi(2));

        fts = fts+ exp(log1pow0 + log1pow1 + rdm1*logts + bico);
        ts = 1. - ts;
    end

    fts = fts*0.5;

    t2 = (ts-0.5) / k;          %   // calculate 1/Ypsilon(ts)
    rgts = exp(-(phideri0 + bico - 0.5 * t2*t2));
    testnum= rh < G_integral * fts * rgts;   % // Bernouilli variate
else
    %   // G > h: can't use sampling in t-domain
     testnum= rh < aux.WalleniusNCHypergeometricpdf(x,n,m,N,odds);
   % testnum= rh < probability(x);
end
end


function bico=lnbico(n,x,N,m)
x2 = n-x; m2 = N-m;
mFac=logfactorial(m)+logfactorial(m2);
% xLastBico =-99;

% natural log of binomial coefficients.
% returns lambda = log(m!*x!/(m-x)!*m2!*x2!/(m2-x2)!)

xFac = logfactorial(x) + logfactorial(x2) + logfactorial(m-x) + logfactorial(m2-x2);
bico=mFac-xFac;
end


function [rd,E,r,w]=findpars(n,x,N,m,omega)
%  // calculate d, E, r, w

xx=[x, n-x];

j=0;
if omega>1
    oo=[1, 1/omega];
else
    oo=[omega, 1];
end
dd=oo(1)*(m-x)+oo(2)*(N-m-xx(2));
d1=1/dd;
E= (oo(1)*m + oo(2)*(N-m)) * d1;

rr=1;
if rr<=d1
    rr=1.2*d1;
    % Newton-Raphson iteration to find r
end

condexit=0;
while condexit==0
    lastr = rr;
    rrc = 1. / rr;
    z = dd - rrc;
    zd = rrc * rrc;
    for i=[1 2]

        rt = rr * oo(i);
        if (rt < 100.)
            % avoid overflow if rt big
            %r21 = 1-2^rt;
            % r2=2^rt;
            [r21,r2]=pow2_1(rt);

            % pow2_1(rt, &r2);         // r2=2^r, r21=1.-2^r
            a = oo(i) / r21;          %      // omegai/(1.-2^r)
            b = xx(i) * a;            %     // x*omegai/(1.-2^r)
            z  = z+ b;
            zd = zd+ b * a * log(2) * r2;
        end
    end

    if zd == 0
        error("FSDA:WalleniusNCHypergeometricrnd:WrgIter",'cannot find r in function CWalleniusNCHypergeometric::findpars');
    end

    rr = rr- z / zd;
    if rr <= d1
        rr = lastr * 0.125 + d1*0.875;
    end
    j=j+1;

    if j== 70
        error("FSDA:WalleniusNCHypergeometricrnd:WrgIter",'convergence problem searching for r in function CWalleniusNCHypergeometric::findpars');
    end


    if abs(rr-lastr) < rr * 1e-6
        condexit=1;
    end
end

if omega > 1
    dd = dd*omega;
    rr = rr*oo(2);
end


r = rr;
rd = rr * dd;


% find peak width

ro = r * omega;
if (ro < 300)                      %  avoid overflow
    [k1,~] = pow2_1(ro);
    k1 = -1. / k1;
    k1 = omega*omega*(k1+k1*k1);
else
    k1 = 0;
end

if r < 300                        %  avoid overflow
    [k2,~]= pow2_1(r);
    k2 = -1. / k2;
    k2 = (k2+k2*k2);
else
    k2 = 0.;
end

phi2d = -4.*r*r*(x*k1 + (n-x)*k2);
if (phi2d >= 0.)
    error("FSDA:WalleniusNCHypergeometricrnd:WrgIter",'peak width undefined in function CWalleniusNCHypergeometric::findpars')
    %        /* wr = r = 0.; */

else
    wr = sqrt(-phi2d);
    w = 1/wr;
end
end
% end of findpars

% Auxiliary functions
function [rt, r2]=pow2_1(r)
r2=2^r;
rt=1-r2;
end

function z=log1mx(x,~)
z=log1p(-x);
% z=log(1-x);
end

function z=log1pow(q, x)
% // calculate log((1-e^q)^x) without loss of precision
% z=log((1-exp(q))^x);

if (abs(q) > 0.1)
    y = exp(q);
    y1 = 1 - y;
else
    %  // expand 1-e^q = -summa(q^n/n!) to avoid loss of precision
    y1 = 0;
    qn = 1; i = 1;    ifac = 1;
    y2=Inf;
    while (y1 ~= y2)
        y2 = y1;
        qn = qn*q;  ifac = ifac*i; i=i+1;
        y1 = y1- qn / ifac;
    end
    y = 1. - y1;
end

if (y > 0.1)
    %  // (1-y)^x calculated without problem
    z= x * log(y1);

else
    % // expand ln(1-y) = -summa(y^n/n)
    y1 =1;  i = 1;  z1 = 0;
    z2=inf;
    while (z1 ~= z2)
        z2 = z1;
        y1 = y1*y;
        z1 = z1- y1 / i; i=i+1;
    end
    z =x * z1;
end
end


%FScategory:ProbDist
