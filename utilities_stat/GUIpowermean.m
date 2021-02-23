function out = GUIpowermean(x, s, w)
%GUIpowermean shows the necessary calculations to obtain the power (generalized) mean in a GUI.
%
%
%<a href="matlab: docsearchFS('GUIpowermean')">Link to the help function</a>
%
%
%  The power mean (also known as generalized mean, or Holder mean or
%  Kolmogorov-Negumo function of the mean), is an abstraction of the
%  Pythagorean means. It includes the harmonic $(s=-1)$, geometric
%  ($s\rightarrow 0$), and arithmetic ($s=1$), quadratic ($s=2$), cubic
%  mean ($s=3$). When $s \rightarrow -\infty$ the power mean tends to
%  $x_{min}$ and when $s\rightarrow +\infty$ the power mean tends to
%  $x_{max}$. http://en.wikipedia.org/wiki/Generalized_mean
%
%  Required input arguments:
%
%     x : vector of numeric data. Vector.
%           Vector containing strictly numerical data.
%           Data Types - double
%
%    s  : power indicator. Scalar.
%         Power indicator for the desired mean ($s=- \inf= x_{min}$; $s=-1$
%         harmonic mean; $s=0$ geometric mean; $s=1$ = arithmetic mean; $s=2$
%         quadratic mean; $s=3$ cubic mean, $s=+\inf$ = $x_{max}$).
%           Example - 1:10
%           Data Types - double
%
% Optional input arguments:
%
%    w  : weights. Vector.
%         Vector of the same length of x containing the weights assigned to
%         each obsearvation. If w is not supplied we assume that all
%         observations have weight equal to 1.
%           Example - 1:10
%           Data Types - double

%
% Output:
%
%    out = detailed output to compute the index. Table. 
%          Table with n+1 rows (where n is the length of x) containing
%          what is shown in the GUI. Last row contains the totals.
%
% See also: GUIvar, GUIstd, GUIquantile, GUIconcentration
%
% References:
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('GUIpowermean')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%

%{
    %% Example of use of the cubic mean.
    x=[60 30 50 50  80 35000  80 95];
    out=GUIpowermean(x,3);
%}

%{
    %% Example of use of weighted geometric mean.
    x=[1.01 1.02 1.03];
    w=[2 6 4];
    out=GUIpowermean(x,0,w);
%}

%% Beginning of code
x=x(:);
n=length(x);
seq=(1:n)';

xs=x.^s;
sstr=num2str(s);

if nargin<3 % unweighted generalized mean
    if abs(s)>1e-7
        
        header={'i' 'x_i' ['x_i^{' sstr '}']};
    else
        header={'i' 'x_i' '\log(x_i)'};
        xs=log(x);
    end
    
    corpus=[seq, x, xs];
    sumdev2=sum(xs);
    
    footer=[NaN sum(x) sumdev2];
    den=n;
    if abs(s)>1e-7
        mx=(sum(xs)/den)^(1/s);
    else
        mx=exp(mean(xs));
    end
    strtitle=['Details of generalized mean calculation when $s$=' sstr];
else % weighted variance
    w=w(:);
    if abs(s)>1e-7
        header={'i' 'x_i' ['x_i^{' sstr '}' ] 'w_i' ['x_i^' sstr 'w_i']};
    else
        header={'i' 'x_i' '\log(x_i)'  'w_i'  ' w_i \log(x_i)'}   ;
        xs=log(x);
    end
    xsw=xs.*w;
    sumw=sum(w);
    dev2=xs.*w;
    corpus=[seq, x, xs, w, xsw];
    sumdev2=sum(dev2);
    if abs(s)>1e-7
        mx=(sum(dev2)/sumw)^(1/s);
    else
        mx=exp(sum(dev2)/sumw);
    end
    footer=[NaN NaN NaN sum(w) sum(xsw) ];
    strtitle=['Details of weighted generalized mean calculation when $s$=' sstr];
    den=sumw;
end


str=strForSchool(header, corpus, footer);

out=array2table([corpus;footer],'VariableNames',header);

fs=14;
dim = [.2 .80 0.1 0.1];
figure('Position',[100 100 1000 600],'Units','normalized');
annotation('textbox',dim,'FitBoxToText','on','String',str,'Interpreter','latex','FontSize',fs);

dim = [.2 .9 0.1 0.1];
fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strtitle,'Interpreter','latex','FontSize',fs1);

dim = [.2 .05 0.1 0.1];

if nargin<3
    if any(isinf(xs))
        strfin=[' \it $\lim_{s \rightarrow + \infty} M_s(X)= \lim_{s \rightarrow + \infty} \left( \frac{\sum_{i=1}^n x_i^s}{n} \right)^{1/s}=x_{\max} =' num2str(max(x)) '$'];
    elseif max(abs(xs))<1e-12 && min(abs(xs))>=0
        strfin=[' \it $\lim_{s \rightarrow - \infty} M_s(X)= \lim_{s \rightarrow - \infty} \left( \frac{\sum_{i=1}^n x_i^s}{n} \right)^{1/s}=x_{\min} =' num2str(min(x)) '$'];
    else
    if abs(s)>1e-7
        strfin=[' \it M$_s(X)= \left( \frac{\sum_{i=1}^n x_i^s}{n} \right)^{1/s}= (\frac{'  num2str(sumdev2) '}{' num2str(den) '})^{1/' sstr '}=' num2str(mx) '$'];
    else
        strfin=[' \it M$_s(X)= \exp \left( \frac{\sum_{i=1}^n \log(x_i)}{n} \right)= \exp (\frac{'  num2str(sumdev2) '}{' num2str(den) '})=' num2str(mx) '$'];
    end
    end
else
    if abs(s)>1e-7
        strfin=[' \it M$_s(X)= \left( \frac{\sum_{i=1}^n x_i^s w_i}{\sum_{i=1}^n w_i} \right)^{(1/s)}= (\frac{'  num2str(sumdev2) '}{' num2str(den) '})^{1/' sstr '}=' num2str(mx) '$'];
    else
        strfin=[' \it M$_s(X)= \exp \left( \frac{\sum_{i=1}^n w_i \log(x_i)}{ \sum_{i=1}^n w_i} \right)= \exp (\frac{'  num2str(sumdev2) '}{' num2str(den) '})=' num2str(mx) '$'];
    end
end

fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strfin,'Interpreter','latex','FontSize',fs1);

end
%FScategory:GUI