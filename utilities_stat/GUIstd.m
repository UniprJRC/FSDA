function out = GUIstd(x,w)
%GUIstd shows the necessary calculations to obtain the standard deviation in a GUI.
%
%
%<a href="matlab: docsearchFS('GUIstd')">Link to the help function</a>
%
%
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
%    w  : weights. Vector or scalar.
%         Vector of the same length of x containing the weights
%         (frequencies) assigned to each observation or scalar which
%         specifies the normalization to use.
%         If w=1 the denominator of the index corresponds to the number of
%         observations.
%         If w=0 (default), the denominator of the index corresponds to the
%         number of observations minus 1,  If w is not supplied we assume
%         that all observations have weight equal to 1 and the index is
%         normalized with length(x)-1.
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
% See also: GUIvar, GUImad, GUIskewness
%
% References:
% Milioli, M.A., Riani, M., Zani, S. (2019), "Introduzione all'analisi dei dati statistici (Quarta edizione ampliata)". [MRZ]
% Cerioli, A., Milioli, M.A., Riani, M. (2016), "Esercizi di statistica (Quinta edizione)". [CMR]
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('GUIstd')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%
%{
    % Standard deviation of the first 6 natural numbers.
    x=1:6;
    GUIstd(x)
%}

%{
    % Calculation of standard deviation (using n as denominator).
    % See page 9 of [CMR]
    x=[427 492 445 444 476 470];
    GUIstd(x,1)
%}

%{
    %% Example of weighted standard deviation. (See page 128 of [MRZ])
    x=[25:50:175 250 400 750];
    w=[59 69 31 25 8 3 1];
    y=GUIstd(x,w);
%}

%% Beginning of code

x=x(:);
lenx=length(x);
seq=(1:lenx)';

if nargin<2 || (nargin==2 && isscalar(w)) % unweighted standard deviation
    nummeanx=sum(x);
    mx=nummeanx/lenx;
    
    header={'i' 'x_i' '(x_i-M)' '(x_i-M)^2'};
    
    dev2=(x-mx).^2;
    corpus=[seq, x, x-mx, dev2];
    sumdev2=sum(dev2);
    
    footer=[NaN sum(x) 0 sumdev2];
    if nargin<2
        w=0;
    end
    if w==0
        den=lenx-1;
    else
        den=lenx;
    end
    strtitle='Details of standard deviation $(\sigma)$ calculation';
    
else % weighted std
    w=w(:);
    xw=x.*w;
    
    den=sum(w);
    nummeanx=sum(xw);
    mx=nummeanx/den;
    lenx=den;
    header={'i' 'x_i' 'w_i' 'x_i w_i' '(x_i-M)w_i' '(x_i-M)^2 w_i'};
    dev2=((x-mx).^2).*w;
    corpus=[seq, x, w, xw, (x-mx).*w, dev2];
    sumdev2=sum(dev2);
    footer=[NaN NaN sum(w) sum(xw) 0 sumdev2];
    strtitle='Details of weighted standard deviation calculation';
    
end
sigmax=sqrt(sumdev2/den);


str=strForSchool(header, corpus, footer);


out=array2table([corpus;footer],'VariableNames',header);

fs=14;
dim = [.2 .80 0.1 0.1];
figure('Position',[100 100 1000 600],'Units','normalized');
annotation('textbox',dim,'FitBoxToText','on','String',str,'Interpreter','latex','FontSize',fs);

dim = [.7 .7 0.1 0.1];
strmean=['\boldmath{$M$}= $\frac{' num2str(nummeanx) '}{' num2str(lenx) '}=' num2str(mx) '$'];
annotation('textbox',dim,'FitBoxToText','on','String',strmean,'Interpreter','latex','FontSize',fs);

dim = [.2 .9 0.1 0.1];
fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strtitle,'Interpreter','latex','FontSize',fs1);

dim = [.2 .05 0.1 0.1];

var=sumdev2/den;
if nargin<2 || isscalar(w)
    if w==0
        strfin=[' \it $\sigma(X)= \sqrt{ \frac{\sum_{i=1}^n (x_i -M)^2}{n-1}}= \sqrt{\frac{'  num2str(sumdev2) '}{' num2str(den) '}}=\sqrt{' num2str(var) '}='  num2str(sigmax) '$'];
    else
        strfin=[' \it $\sigma(X)= \sqrt{\frac{\sum_{i=1}^n (x_i -M)^2}{n}}= \sqrt{\frac{'  num2str(sumdev2) '}{' num2str(den) '}}=\sqrt{' num2str(var) '}='  num2str(sigmax) '$'];
    end
else
    strfin=[' \it $\sigma(X)= \sqrt{\frac{\sum_{i=1}^n (x_i -M)^2 w_i}{\sum_{i=1}^n w_i}}= \sqrt{ \frac{'  num2str(sumdev2) '}{' num2str(den) '}}=\sqrt{' num2str(var)  '}='  num2str(sigmax) '$'];
end

fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strfin,'Interpreter','latex','FontSize',fs1);

end
%FScategory:GUI