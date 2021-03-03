function out = GUIskewness(x, flag, w)
%GUIskewness shows the necessary calculations to obtain the variance in a GUI.
%
%
%<a href="matlab: docsearchFS('GUIskewness')">Link to the help function</a>
%
%
%  Required input arguments:
%
%     x : vector of numeric data. Vector.
%           Vector containing strictly numerical data.
%           Data Types - double
%
%
%  Optional input arguments:
%
%    flag  : flag for unbiasness. Scalar.
%         If flag=0 the unbiased index is computed.
%         If flag=1 (default) the biased index is computed,
%         meaning that it tends to differ from the population skewness by a
%         systematic amount based on the sample size.
%           Example - 1
%           Data Types - double
%
%    w  : weights. Vector.
%         Vector of the same length of x containing the weights assigned to
%         each observation.
%           Example - 1:10
%           Data Types - double
%
% Output:
%
%    out = detailed output to compute the index. Table. 
%          Table with n+1 rows (where n is the length of x) containing
%          what is shown in the GUI. Last row contains the totals.
%
%
%
%
% See also: GUIstd, GUImad, GUIquantile
%
% References:
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('GUIskewness')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%
%{
    %% Calculation of biased index of skewness.
    % Vector x contains the expenditure of 20 customers in a supermarket.
    x=[20, 110,  31, 125,  40, 128,  50, 140,  65, 150,  72, 175,...
         85, 190,  100, 205,  100, 230,  105 310];
    GUIskewness(x)
%}

%{
    % Calculation of unbiased index of skewness.
    % Vector x contains the expenditure of 20 customers in a supermarket.
    x=[20, 110,  31, 125,  40, 128,  50, 140,  65, 150,  72, 175,...
         85, 190,  100, 205,  100, 230,  105 310];
    GUIskewness(x,0)
%}

%{
    % Example of weighted unbiased index of skewness.
    % Italian grades obtained in a particular university exam by 1000 students.
    x=20:29;
    % w= vector of frequencies
    w=[61    82   102   164   122   122    61   143    82    61];
    GUIskewness(x,0,w)
%}

%{
    % Example of weighted biased index of skewness.
    % Italian grades obtained in a particular university exam by 1000 students.
    x=20:29;
    % w= vector of frequencies
    w=[61    82   102   164   122   122    61   143    82    61];
    GUIskewness(x,1,w)
%}

%% Beginning of code

x=x(:);
n=length(x);
seq=(1:n)';


if nargin<3  % unweighted skewness
    if nargin<2
        flag=1;
    end
    nummeanx=sum(x);
    mx=nummeanx/n;
    
    header={'i' 'x_i' '(x_i-M)' '(x_i-M)^2'  '(x_i-M)^3'};
    
    dev1=(x-mx);
    dev2=dev1.^2;
    dev3=dev1.^3;
    corpus=[seq, x, dev1, dev2, dev3];
    sumdev2=sum(dev2);
    sumdev3=sum(dev3);
    
    footer=[NaN sum(x) 0 sumdev2 sumdev3];
    
    if flag==0
        corr=sqrt(n*(n-1))/(n-2);
        strtitle='Details of skewness calculation (unbiased version) $\gamma_u(X)$';
    elseif flag==1
        corr=1;
        strtitle='Details of skewness calculation (biased version)';
    else
        error(message('FSDA:GUIskewnwss:BadFlagReduction 2nd input arg must be 0 or 1'));
    end
    dim = [.2 .80 0.1 0.1];
    
else
    
    
    % weighted skewness
    w=w(:);
    xw=x.*w;
    sumw=sum(w);
    n=sumw;
    nummeanx=sum(xw);
    mx=nummeanx/sumw;
    header={'i' 'x_i' 'n_i' 'x_i n_i' '(x_i-M)n_i' '(x_i-M)^2 n_i'  '(x_i-M)^3 n_i'};
    dev2=((x-mx).^2).*w;
    dev3=((x-mx).^3).*w;
    corpus=[seq, x, w, xw, (x-mx).*w, dev2 dev3];
    sumdev2=sum(dev2);
    sumdev3=sum(dev3);
    footer=[NaN NaN sum(w) sum(xw) 0 sumdev2 sumdev3];
    
    if flag==0
        corr=sqrt(n*(n-1))/(n-2);
        strtitle='Details of weighted skewness calculation (unbiased version) $\gamma_u(X)$';
    elseif flag==1
        corr=1;
        strtitle='Details of weighted skewness calculation (biased version)';
    else
        error(message('FSDA:GUIskewnwss:BadFlagReduction 2nd input arg must be 0 or 1'));
    end
    dim = [.02 .80 0.1 0.1];
end


str=strForSchool(header, corpus, footer);


out=array2table([corpus;footer],'VariableNames',header);

fs=14;
figure('Position',[100 100 1000 600],'Units','normalized');
annotation('textbox',dim,'FitBoxToText','on','String',str,'Interpreter','latex','FontSize',fs);

dim = [.7 .7 0.1 0.1];
strmean=['\boldmath{$M$}= $\frac{' num2str(nummeanx) '}{' num2str(n) '}=' num2str(mx) '$'];
annotation('textbox',dim,'FitBoxToText','on','String',strmean,'Interpreter','latex','FontSize',fs);

dim = [.2 .9 0.1 0.1];
fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strtitle,'Interpreter','latex','FontSize',fs1);

Fisherind=corr*sqrt(n)*sumdev3/sumdev2^1.5;
nstr=num2str(n);

if nargin<3 % all weights equal to 1
    if flag==0 % unbiased inded
        gammaind=' \frac{\sqrt{n(n-1)}}{n-2} \gamma(X)=  \frac{\sqrt{n(n-1)}}{n-2}  \frac{\sqrt n \sum_{i=1}^n (x_i -M(X))^3 }{\left[    \sum_{i=1}^n (x_i -M(X))^2 \right]^{3/2}}';
        gammaind=['\gamma_{u}(X)= ' gammaind];
        strfin=[' \it $' gammaind '=\frac{\sqrt{' nstr ' \times ' num2str(n-1) '}}{' num2str(n-2) '}  \frac{ \sqrt{' nstr '} \times'  num2str(sumdev3) '}{' num2str(sumdev2) '^{1.5}}=' num2str(Fisherind) '$'];
        dim = [.01 .05 0.1 0.1];
    else  % biased index
        gammaind='\gamma(X)= \frac{\sqrt n \sum_{i=1}^n (x_i -M(X))^3 }{\left[    \sum_{i=1}^n (x_i -M(X))^2 \right]^{3/2}}';
        strfin=[' \it $' gammaind '= \frac{ \sqrt{' nstr '} \times'  num2str(sumdev3) '}{' num2str(sumdev2) '^{1.5}}=' num2str(Fisherind) '$'];
        dim = [.2 .05 0.1 0.1];
    end
else % frequency distribution
    if flag==0 % unbiased index
        gammaind=' \frac{\sqrt{n(n-1)}}{n-2} \gamma(X)=  \frac{\sqrt{n(n-1)}}{n-2}  \frac{\sqrt n \sum_{i=1}^r (x_i -M(X))^3 n_i }{\left[    \sum_{i=1}^r (x_i -M(X))^2 n_i \right]^{3/2}}';
        gammaind=['\gamma_{u}(X)= ' gammaind];
        strfin=[' \it $' gammaind '=\frac{\sqrt{' nstr ' \times ' num2str(n-1) '}}{' num2str(n-2) '}  \frac{ \sqrt{' nstr '} \times'  num2str(sumdev3) '}{' num2str(sumdev2) '^{1.5}}=' num2str(Fisherind) '$'];
        dim = [.01 .05 0.1 0.1];
    else  % biased index
        gammaind='\gamma(X)= \frac{\sqrt n \sum_{i=1}^r (x_i -M(X))^3 n_i}{\left[    \sum_{i=1}^r (x_i -M(X))^2 n_i \right]^{3/2}}';
        strfin=[' \it $' gammaind '= \frac{ \sqrt{' nstr '} \times'  num2str(sumdev3) '}{' num2str(sumdev2) '^{1.5}}=' num2str(Fisherind) '$'];
        dim = [.2 .05 0.1 0.1];
    end
end

fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strfin,'Interpreter','latex','FontSize',fs1);

end
%FScategory:GUI