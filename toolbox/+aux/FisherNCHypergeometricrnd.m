function r = FisherNCHypergeometricrnd(n,m,N,odds, accuracy)
%FisherNCHypergeometricrnd returns a random number from non-central hypergeometric probability density function
if nargin<4
    odds=1;
    accuracy=1e-10;
elseif  nargin<5
    accuracy=1e-10;
else

end


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
         error('FSDA:FisherNCHypergeometricrnd:Wrginpt',"Not enough items with nonzero weight in function WalleniusNCHyp");
    end
end


if (odds == 1)
    % // use hypergeometric function if odds == 1
    % M. K. N
    r = hygernd(N, m, m);
    return
end

% symmetry transformations
fak = 1;  addd = 0;
if (m > N/2)
    % // invert m
    m = N - m;
    fak = -1;  addd = n;
end

if (n > N/2)
    % // invert n
    n = N - n;
    addd = addd+ fak * m;
    fak = - fak;
end


if (n > m)
    % // swap n and m
    x = n;  n = m;  m = x;
end

% // cases with only one possible result end here
if (n == 0 || odds == 0.)
    r=addd;
    return
end

if (fak == -1)
    % // reciprocal odds if inverting
    odds = 1. / odds;
end

% // choose method
if (n < 30 && N < 1024 && odds > 1.E-5 && odds < 1.E5)
    % // use inversion by chop down method
    x = FishersNCHypInversion(n, m, N, odds);
else
    % // use ratio-of-uniforms method
    x = FishersNCHypRatioOfUnifoms(n, m, N, odds);
end

% // undo symmetry transformations
r = x * fak + addd;
end

function x=FishersNCHypInversion(n, m, N, odds)
%  Subfunction for FishersNCHyp distribution.
%  Implements Fisher's noncentral hypergeometric distribution by inversion
%  method, using chop-down search starting at zero.

% Valid only for 0 <= n <= m <= N/2.
% Without overflow check the parameters must be limited to n < 30, N < 1024,
% and 1.E-5 < odds < 1.E5. This limitation is acceptable because this method
% is slow for higher n.

%The execution time of this function grows with n.

%See the file nchyp.pdf for theoretical explanation.
% int32_t x;                          // x value
% double f;                           // scaled function value
% double sum;                         // scaled sum of function values
% double a1, a2, b1, b2, f1, f2;      // factors in recursive calculation
% double u;                           // uniform random variate

L = N - m - n;  % derived parameter


% // f(0) is set to an arbitrary value because it cancels out.
% // A low value is chosen to avoid overflow.
fnc_f0 = 1.E-100;

% // calculate summation of e(x), using the formula:
% // f(x) = f(x-1) * (m-x+1)*(n-x+1)*odds / (x*(L+x))
% // All divisions are avoided by scaling the parameters
sumd =  fnc_f0;
f =  fnc_f0;
fnc_scale = 1.;
a1 = m;  a2 = n;  b1 = 1;  b2 = L + 1;
for x = 1:n
    f1 = a1 * a2 * odds;
    f2 = b1 * b2;
    a1=a1-1;
    a2=a2-1;
    b1=b1+1;
    b2=b2+1;
    f = f*f1;
    sumd = sumd* f2;
    fnc_scale = fnc_scale * f2;
    sumd = sumd+ f;
end
fnc_f0 = fnc_f0* fnc_scale;
fnc_scale = sumd;


% // uniform random
u = unifrnd(0,1) * fnc_scale;

% // recursive calculation:
% // f(x) = f(x-1) * (m-x+1)*(n-x+1)*odds / (x*(L+x))
f = fnc_f0;  x = 0;  a1 = m;  a2 = n;  b1 = 0;  b2 = L;
while x<n
    u = u-f;
    if (u <= 0)
        break
    end

    x=x+1;  b1=b1+1;  b2=b2+1;
    f = f*a1 * a2 * odds;
    u = u*b1 * b2;
    a1=a1-1;  a2=a2-1;
end

end


function k = FishersNCHypRatioOfUnifoms(n, m, N, odds)
 %{ 
   Subfunction for FishersNCHyp distribution. 
   Valid for 0 <= n <= m <= N/2, odds != 1

   Fisher's noncentral hypergeometric distribution by ratio-of-uniforms 
   rejection method.

   The execution time of this function is almost independent of the parameters.
 %}
   % int32_t L;                          
   % int32_t mode;                       // mode
   % double mean;                        // mean
   % double variance;                    // variance
   % double x;                           // real sample
   % int32_t k;                          // integer sample
   % double u;                           // uniform random
   % double lf;                          // ln(f(x))
   % double AA, BB, g1, g2;              // temporary

   L = N - m - n; % // N-m-n

   
      % // find approximate mean
      AA = (m+n)*odds+L; BB = sqrt(AA*AA - 4*odds*(odds-1)*m*n);
      mea = (AA-BB)/(2*(odds-1));

      % // find approximate variance
      AA = mea * (m-mea); BB = (n-mea)*(mea+L);
      vari = N*AA*BB/((N-1)*(m*BB+(n+L)*AA));

      % // compute log(odds)
      fnc_logb = log(odds);

      % // find center and width of hat function
      fnc_a = mea + 0.5;
      fnc_h = 1.028 + 1.717*sqrt(vari+0.5) + 0.032*abs(fnc_logb);

      % find safety bound
      fnc_bound = floor(mea + 4.0 * fnc_h);
      if fnc_bound > n 
          fnc_bound = n;
      end

      % // find mode
      mode = floor(mea);
      g1 =(m-mode)*(n-mode)*odds;
      g2 =(mode+1)*(L+mode+1);
      if (g1 > g2 && mode < n) 
          mode=mode+1;
      end

      % // value at mode to scale with:
      fnc_lfm = mode * fnc_logb - fc_lnpk(mode, L, m, n);
   

   while 1>0
      u = unifrnd(0,1);
      if (u == 0) 
          continue                   %  // avoid divide by 0
      end

      u2=unifrnd(0,1);

      x = fnc_a + fnc_h * (u2-0.5)/u;
      if (x < 0. || x > 2E9) 
          continue;         %  // reject, avoid overflow
      end

      k = floor(x);                         %  // truncate
      if (k > fnc_bound) 
          continue;              %  // reject if outside safety bound
      end

      lf = k*fnc_logb - fc_lnpk(k,L,m,n) - fnc_lfm; % // compute function value
      if (u * (4.0 - u) - 3.0 <= lf) 
          break     % // lower squeeze accept
      end

      if (u * (u-lf) > 1.0) 
          continue;            % // upper squeeze reject
      end

      if (2.0 * log(u) <= lf) 
          break             %v // final acceptance
      end
   end
end


function lnpk=fc_lnpk(k, L, m, n) 
   % // subfunction used by hypergeometric and Fisher's noncentral hypergeometric distribution
   lnpk=logfactorial(k) + logfactorial(m - k) + logfactorial(n - k) + logfactorial(L + k);
   end

%FScategory:ProbDist