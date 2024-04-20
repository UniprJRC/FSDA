
function Wpdf = CMultiWalleniusNCHypergeometricpdf(x,m,omega, accuracy)
%CMultiWalleniusNCHypergeometricpdf returns multivariate Wallenius' non-central hypergeometric probability density function
% n = number of balls taken
%  m = number of balls for each color
%  total number of balls in the urn = sum(m)
% x = number of successes for each ball
% omega = vector of odds

if nargin<4
    accuracy =1e-09;
end

m=m(:);
omega=omega(:);


colors=length(m);
n=sum(x);
N=sum(m);

if (colors < 3)
    if (colors <= 0)
        Wpdf=1;
        return
    end

    if (colors == 1)
        Wpdf=x(1) == m(1);
        return
    end

    if colors == 2

        if (omega(2) == 0.)
            Wpdf=x(1) == m(1);
            return
        end
       
        odds=omega(1)/omega(2);
        %  x,n,m,N,omega
        Wpdf= aux.WalleniusNCHypergeometricpdf(x(1),n,m(1),N,odds,accuracy);
        return
    end
end


j=0;
em=0;

central = 1;
for i  = 0:1:colors-1
    if x(i+1) > m(i+1) || x(i+1) < 0 || x(i+1) < n - N + m(i+1)
        Wpdf=0;
        return
    end

    if (x(i+1) > 0)
        j=j+1;
    end

    if (omega(i+1) == 0. && x(i+1)~=0)
        Wpdf=0;
        return
    end

    if (x(i+1) == m(i+1) || omega(i+1) == 0.)
        em=em+1;
    end

    if i > 1 && omega(i+1) ~= omega(i)
        central = 0;
    end

end

if (n == 0 || em == colors)
    Wpdf=1;
    return
end

if central ==1

    %// All omega's are equal.
    % // This is multivariate central hypergeometric distribution
    p = 1;
    sx=n; sm=N;
    for i = 1:colors
        % // Use univariate hypergeometric (usedcolors-1) times
        p = p* hygepdf(x(i), sm, m(i),sx);
        sx = sx- x(i);
        sm = sm -m(i);
    end
    Wpdf=p;
    return
end

if j == 1
    Wpdf=binoexpand(x,m,omega);
    return
end


[rd,~,r,w]=findpars(x,m,omega);
%  // calculate rd, E, r, w
Wpdf=integrate(w,rd,r, x,m,omega );

end

function [rd,E,r,w]= findpars(x,m,omega)

r=1;
%// calculate r, w, E
%  // calculate d, E, r, w
colors=length(omega);

% // find highest omega
omax=max(omega);
omaxr = 1. / omax;
dd =0;
E = 0;
omeg=omega;
for i = 1:colors
    % // scale weights to make max = 1
    omeg(i) = omega(i) * omaxr;
    % // calculate d and E
    dd = dd+omeg(i) * (m(i)-x(i));
    E  = E+ omeg(i) * m(i);
end
dr = 1. / dd;
E = E*dr;
rr = r * omax;
if (rr <= dr)
    rr = 1.2 * dr;  % // initial guess
end

% // Newton-Raphson iteration to find r
lastr=Inf;
j=1;
while abs(rr-lastr) > rr * 1E-5
    lastr = rr;
    rrc = 1. / rr;
    z = dd - rrc;                   % // z(r)
    zd = rrc * rrc;                  % // z'(r)
    for i=1:colors
        rt = rr * omeg(i);
        if (rt < 100. && rt > 0.)    %// avoid overflow and division by 0
            [r21,r2] = pow2_1(rt);    % // r2=2^r, r21=1.-2^r
            a = omeg(i) / r21;        % // omegai/(1.-2^r)
            b = x(i) * a;             % // x*omegai/(1.-2^r)
            z  = z+b;
            zd =zd+ b * a * r2 * log(2);
        end
    end
    if (zd == 0)
        error('FSDA:CMultiWalleniusNCHypergeometricpdf:WrgIter',"can't find r in function CMultiWalleniusNCHypergeometric::findpars");
    end

    rr = rr- z / zd;                   %  // next r
    if (rr <= dr)
        rr = lastr * 0.125 + dr * 0.875;
    end
    j=j+1;
    if j == 70
        error("convergence problem searching for r in function CMultiWalleniusNCHypergeometric::findpars");
    end
end
rd = rr * dd;
r = rr * omaxr;

% // find peak width
phi2d = 0.;
for i=1:colors
    ro = rr * omeg(i);
    if (ro < 300 && ro > 0.)        % // avoid overflow and division by 0
        k1 = pow2_1(ro);
        k1 = -1. / k1;
        k1 = omeg(i) * omeg(i) * (k1 + k1*k1);
    else
        k1 = 0;
    end
    phi2d = phi2d + x(i) * k1;
end
phi2d = -phi2d *(4 * rr * rr);
if (phi2d > 0.)
    error('FSDA:CMultiWalleniusNCHypergeometricrnd:WrgInpt',"peak width undefined in function CMultiWalleniusNCHypergeometric::findpars");
end
wr = sqrt(-phi2d);
w = 1. / wr;
end

function [rt, r2]=pow2_1(r)
r2=2^r;
rt=1-r2;
end


%% Beginning of integrate
function  Wpdf=integrate(w,rd,r,  x,m,omega, accuracy)

if nargin<7
    accuracy=1e-10;
end

bico=lnbico(x,m,omega);

if w < 0.02
    if accuracy < 1E-9
        s1=0.5;
    else
        s1=1;
    end

    delta = s1 * w;                      % // integration steplength
    ta = 0.5 + 0.5 * delta;
    sumd = integrate_step(1.-ta, ta,rd,r,omega,x,m,bico);      % // first integration step around center peak
    tb=0;
    while (tb<1)
        tb = ta + delta;
        if (tb > 1.)
            tb = 1;
        end

        s  = integrate_step(ta, tb,rd,r,omega,x,m,bico);       % // integration step to the right of peak
        s = s+integrate_step(1.-tb,1.-ta,rd,r,omega,x,m,bico);  %// integration step to the left of peak
        sumd = sumd+s;
        if (s < accuracy * sumd)
            % stop before interval finished if accuracy reached
            break
        end

        ta = tb;
        if (tb > 0.5 + w)
            delta = delta*2;    %  // increase step length far from peak
        end
    end
else
    % difficult situation. Step length determined by inflection points
    sumd=0;
    for t1=[0 0.5]

        t2=t1+0.5;

        tinf = search_inflect(t1, t2, x,m,omega, rd,r);     % find inflection point

        delta = tinf - t1;
        if (delta > t2 - tinf)
            delta = t2 - tinf; % // distance to nearest endpoint
        end

        delta = delta/7;

        if (delta < 1E-4)
            delta = 1E-4;
        end

        delta1 = delta;
        % // integrate from tinf forwards to t2
        ta = tinf;
        tb=0;

        while (tb < t2)
            tb = ta + delta1;
            if (tb > t2 - 0.25*delta1)
                tb = t2; % // last step of this subinterval
            end

            s = integrate_step(ta, tb,rd,r,omega,x,m,bico);        %  // integration step
            sumd = sumd+s;
            delta1 = delta1*2;                 %        // double steplength
            if (s < sumd * 1E-4)
                delta1 = delta1* 8;   %  large step when s small
            end
            ta = tb;

        end

        if tinf
            % // integrate from tinf backwards to t1
            tb = tinf;
            while (ta > t1)
                ta = tb - delta;
                if (ta < t1 + 0.25*delta)
                    ta = t1; % // last step of this subinterval
                end
                s = integrate_step(ta, tb,rd,r,omega,x,m,bico);       % // integration step
                sumd = sumd+ s;
                delta = delta * 2;                       % // double steplength
                if (s < sumd * 1E-4)
                    delta = delta* 8;  % // large step when s small
                end
                tb = ta;
            end
        end
    end
end
Wpdf= sumd*rd;
end


function bico=lnbico(x,m,omega)
% // natural log of binomial coefficients
bico = 0;
colors=length(m);
for i=1:colors
    if x(i) < m(i) && omega(i)~=0
        bico = bico + logfactorial(m(i)) - logfactorial(x(i)) - logfactorial(m(i)-x(i));
    end
end
end


function t=search_inflect(t_from, t_to, x,m,omega, rd,r)
xx=x;
colors=length(m);
rdm1 = rd - 1.;
if (t_from == 0 && rdm1 <= 1.)
    t=0;
    return % //no inflection point
end

t = 0.5 * (t_from + t_to);
zeta=zeros(colors,3,3);
rho=zeros(colors,1);
for i=1:colors
    rho(i)=r*omega(i);
    % // calculate zeta coefficients
    zeta(i,1,1) = rho(i);
    zeta(i,1,2) = rho(i) * (rho(i) - 1.);
    zeta(i,2,2) = rho(i) * rho(i);
    zeta(i,1,3) = zeta(i,1,2) * (rho(i) - 2.);
    zeta(i,2,3) = zeta(i,1,2) * rho(i) * 3.;
    zeta(i,3,3) = zeta(i,2,2) * rho(i) * 2.;
end


iter = 0;
sele=repmat([0;0;1;1],20,1);

t1=Inf;
while (abs(t - t1) > 1e-5)
    t1 = t;

    tr = 1. / t;
    log2t = log(t)*(1./log(2));
    % phi(1) = phi(2] = phi[3] = 0.;
    phi=zeros(3,1);
    for i=1:colors            % // ca lculate first 3 derivatives of phi(t)
        [q1,q] = pow2_1(rho(i)*log2t);
        q =q/ q1;
        phi(1) = phi(1) -xx(i) * zeta(i,1,1) * q;
        phi(2) = phi(2) - xx(i) * q * (zeta(i,1,2) + q * zeta(i,2,2));
        phi(3) = phi(3) - xx(i) * q * (zeta(i,1,3) + q * (zeta(i,2,3) + q * zeta(i,3,3)));
    end

    phi(1) = phi(1)+rdm1;
    phi(2) = phi(2) - rdm1;
    phi(3) = phi(3) + 2. * rdm1;
    phi(1) = phi(1) *tr;
    phi(2) = phi(2)*tr * tr;
    phi(3) = phi(3)* tr * tr * tr;
    % method = (iter & 2) >> 1;       %  // alternate between the two methods ALDO
    method=sele(iter+1);

    Z2 = phi(1)*phi(1) + phi(2);
    Zd = method*phi(1)*phi(1)*phi(1) + (2.+method)*phi(1)*phi(2) + phi(3);

    if (t < 0.5)
        if (Z2 > 0)
            t_from = t;
        else
            t_to = t;
        end
        if (Zd >= 0)
            % // use binary search if Newton-Raphson iteration makes problems
            % t = (t_from ? 0.5 : 0.2) * (t_from + t_to);    % ALDO
            t= t_from;
        else  % Newton-Raphson iteration
            t = t- Z2 / Zd;
        end

    else
        if (Z2 < 0)
            t_from = t;
        else
            t_to = t;
        end

        if (Zd <= 0)
            % // use binary search if Newton-Raphson iteration makes problems
            t = 0.5 * (t_from + t_to);
        else
            % // Newton-Raphson iteration
            t = t - Z2 / Zd;
        end

    end
    if (t >= t_to)
        t = (t1 + t_to) * 0.5;
    end

    if (t <= t_from)
        t = (t1 + t_from) * 0.5;
    end

    iter=iter+1;

    if iter > 20
        error('FSDA:CMultiWalleniusNCHypergeometricrnd:WrgIter',"Search for inflection point failed in function CWalleniusNCHypergeometric::search_inflect");
    end
end
end


function  s = integrate_step(ta, tb,rd,r,omega,x,m,bico)
IPOINTS=8;
colors=length(m);
xval=[-0.960289856498,-0.796666477414,-0.525532409916,-0.183434642496,0.183434642496,0.525532409916,0.796666477414,0.960289856498];
weights= [0.10122853629,0.222381034453,0.313706645878,0.362683783378,0.362683783378,0.313706645878,0.222381034453,0.10122853629];

delta = 0.5 * (tb - ta);
ab = 0.5 * (ta + tb);
rdm1 = rd - 1.;
sumd = 0;

for j=1:IPOINTS
    tau = ab + delta * xval(j);
    ltau = log(tau);
    taur = r * ltau;
    y=0;
    for i = 1:colors
        if omega(i)~=0
            % // possible loss of precision due to subtraction here:
            y = y+log1pow(taur*omega(i),x(i));
        end
    end
    y=y+ rdm1*ltau + bico;
    if (y > -50)
        sumd = sumd+ weights(j) * exp(y);
    end
end
s=delta * sumd;
end


function z=log1mx(x,~)
z=log(1-x);
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


function wpdf=binoexpand(x,m,omega)
   %// binomial expansion of integrand
   % // only implemented for x(i) = 0 for all but one i
   W=sum(omega.*m);
   n=sum(x);
   % total weight
   j=find(x>0,1);
   % if (k > 1) FatalError("More than one x[i] nonzero in CMultiWalleniusNCHypergeometric::binoexpand");
   fal1=FallingFactorial(m(j),n);
   fal2=FallingFactorial(W/omega(j),n);
   wpdf=exp(fal1- fal2);
end


function z=FallingFactorial(a,b)
% // calculates ln(a*(a-1)*(a-2)* ... * (a-b+1))

if (b < 30 && round(b) == b && a < 1E10)
    %  direct calculation
    % double f = 1.;
    %for (int i = 0; i < b; i++) f *= a--;
    %return log(f);
    z=log(prod(a:-1:(a-b+1)));

elseif (a > 100*b && b > 1)

    % % // combine Stirling formulas for a and (a-b) to avoid loss of precision
    ar = 1./a;
    cr = 1./(a-b);
    % % // calculate -log(1-b/a) by Taylor expansion
    s = 0;
    n = 1;
    ba = b*ar;
    f = ba;
    lasts=Inf;
    while (s ~= lasts)
        lasts = s;
        s = s+f/n;
        f = f*ba;
        n=n+1;
    end

    z= (a+0.5)*s + b*log(a-b) - b + (1./12.)*(ar-cr);
    %    //- (1./360.)*(ar*ar*ar-cr*cr*cr)

else
    % // use LnFacr function
    lna=LnFacr(a);
    lnab=LnFacr(a-b);
    z=lna - lnab;
end
end



function f=LnFacr(x)
% log factorial of non-integer x
ix = round(x);
if (x == ix)
    f=logfactorial(ix);
    % // x is integer
else
    D=1;
    C0 =  0.918938533204672722;      % // ln(sqrt(2*pi))
    C1 =  1/12;
    C3 = -1/360;
    C5 =  1/1260;
    C7 = -1/1680;
    if (x < 6)
        if (x == 0 || x == 1)
            f=0;
            return
        else
            while (x < 6)
                D = D*x;
            end
        end
    end
    r  = 1. / x;  r2 = r*r;
    f = (x + 0.5)*log(x) - x + C0 + r*(C1 + r2*(C3 + r2*(C5 + r2*C7)));
    if D~=1
        f = f-  log(D);
    end
end
end
