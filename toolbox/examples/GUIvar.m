function out = GUIvar(x,w)
%GUIvar shows the necessary calculations to obtain the variance in a GUI.
%
%
%<a href="matlab: docsearchFS('GUIvar')">Link to the help function</a>
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
%         Vector of the same length of x containing the weights assigned to
%         each observation or scalar which specifies the normalization to use
%         If w=1 the denominator of the index corresponds to the number of
%         observations.
%         If w=0 (default), the denominator of the index corresponds to the
%         number of observations minus 1. If w is not supplied we assume that
%         all observations have weight equal to 1 and the index is normalized with length(x)-1.
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
% Milioli, M.A., Riani, M., Zani, S. (2019), "Introduzione all'analisi dei dati statistici (Quarta edizione ampliata)". [MRZ]
%
% Cerioli, A., Milioli, M.A., Riani, M. (2016), "Esercizi di statistica (Quinta edizione)". [CMR]
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('GUIvar')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%
%{
    % Variance of the first 6 natural numbers.
    x=1:6;
    out=GUIvar(x);
%}

%{  
    % Calculation of variance (using n as denominator).
    % See page 9 of [CMR]
    x=[427 492 445 444 476 470];
    y=GUIvar(x,1);
%}

%{
    %% Example of weighted variance.
    % See page 128 of [MRZ]
    x=[25:50:175 250 400 750];
    w=[59 69 31 25 8 3 1];
    out=GUIvar(x,w);
%}

%% Beginning of code

x=x(:);
lenx=length(x);
seq=(1:lenx)';

if nargin<2 || (nargin==2 && isscalar(w)) % unweighted variance
    nummeanx=sum(x);
    mx=nummeanx/lenx;
    
    header={'i' 'x_i' '(x_i-M)' '(x_i-M)^2'};
    
    dev2=(x-mx).^2;
    corpus=[seq, x, x-mx, dev2];
    sumdev2=sum(dev2);%     if length(x)>17
    %         corpus=corpus([1:8 lenx-7:lenx],:);
    %     end
    
    footer=[NaN sum(x) 0 sumdev2];
    if nargin<2
        w=0;
    end
    if w==0
        den=lenx-1;
    else
        den=lenx;
    end
    strtitle='Details of variance calculation';
else % weighted variance
    w=w(:);
    xw=x.*w;
    sumw=sum(w);
    lenx=sumw;
    nummeanx=sum(xw);
    mx=nummeanx/sumw;
    header={'i' 'x_i' 'w_i' 'x_i w_i' '(x_i-M)w_i' '(x_i-M)^2 w_i'};
    dev2=((x-mx).^2).*w;
    corpus=[seq, x, w, xw, (x-mx).*w, dev2];
    sumdev2=sum(dev2);
    footer=[NaN NaN sum(w) sum(xw) 0 sumdev2];
    strtitle='Details of weighted variance calculation';
    den=sumw;
end


str=strForSchool(header, corpus, footer);




fs=14;
dim = [.2 .80 0.1 0.1];
figure('Position',[100 100 1000 600],'Units','normalized');
% Make sure that that figure is also visible inside .mlx files
scatter([],[]);
axis('off')
set(gcf,'Visible','on')
annotation('textbox',dim,'FitBoxToText','on','String',str,'Interpreter','latex','FontSize',fs);

dim = [.7 .7 0.1 0.1];
 if verLessThanFS('25.1')
     strmean=['\boldmath{$M$}= $\frac{' num2str(nummeanx) '}{' num2str(lenx) '}=' num2str(mx) '$'];
 else
     strmean=['$M = \frac{' num2str(nummeanx) '}{' num2str(lenx) '}=' num2str(mx) '$'];
 end

annotation('textbox',dim,'FitBoxToText','on','String',strmean,'Interpreter','latex','FontSize',fs);

dim = [.2 .9 0.1 0.1];
fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strtitle,'Interpreter','latex','FontSize',fs1);

dim = [.2 .05 0.1 0.1];

var=sumdev2/den;
if nargin<2 || isscalar(w)
    if w==0
        strfin=[' \it VAR$(X)= \frac{\sum_{i=1}^n (x_i -M)^2}{n-1}= \frac{'  num2str(sumdev2) '}{' num2str(den) '}=' num2str(var) '$'];
    else
        strfin=[' \it VAR$(X)= \frac{\sum_{i=1}^n (x_i -M)^2}{n}= \frac{'  num2str(sumdev2) '}{' num2str(den) '}=' num2str(var) '$'];
    end
else
    strfin=[' \it VAR$(X)= \frac{\sum_{i=1}^n (x_i -M)^2 w_i}{\sum_{i=1}^n w_i}= \frac{'  num2str(sumdev2) '}{' num2str(den) '}=' num2str(var) '$'];
end

fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strfin,'Interpreter','latex','FontSize',fs1);

out=struct;
out.data=array2table([corpus;footer],'VariableNames',header);
out.var=var;
end
%FScategory:GUI