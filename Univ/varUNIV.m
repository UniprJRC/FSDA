function out = varUNIV(x,w)
%varUNIV shows the necessary calculations to obtain the variance in a GUI.
%
%
%<a href="matlab: docsearchFS('varUNIV')">Link to the help function</a>
%
%
%
%
%  Required input arguments:
%
%     x : vector of numeric data. Vector.
%           Vector containg strictly numrical data
%           Data Types - double
% 
%
%  Optional input arguments:
%
%    w  : weights. Vector.
%         Vector of the same legth of x containing the weights assignied to each obsearvation.
%         If w is not supplied we assume that all observations have weight euqual to 1.
%           Example - 1:10
%           Data Types - double
%
% Output:
%
%    out = detailed output to compute the index. Table. Table with n+1 rows (where n is the length of x) containing the
%          what is shown in the GUI. Last row contains the totals.
%
%
%
%
% See also: MixSim, restreigen, restrdeter
%
% References:
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('varUNIV')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%
%{
    % Variance of the first 6 natural numbers.
    x=1:6;
    varUNIV(x)
%}

%{
    % Example of weighted variance.
    x=[25:50:175 250 400 750];
    w=[59 69 31 25 8 3 1];
    y=varUNIV(x,w);
%}

%% Beginning of code

x=x(:);
lenx=length(x);
seq=(1:lenx)';

if nargin<2 % unweighted variance
    mx=mean(x);
    header={'i' 'x_i' '(x_i-M)' '(x_i-M)^2'};
    
    dev2=(x-mx).^2;
    corpus=[seq, x, x-mx, dev2];
    sumdev2=sum(dev2);%     if length(x)>17
%         corpus=corpus([1:8 lenx-7:lenx],:);
%     end        

    footer=[NaN sum(x) 0 sumdev2];
    den=lenx-1;
    strtitle='Details of variance calculation';
else % weighted variance
    w=w(:);
    xw=x.*w;
    sumw=sum(w);
    mx=sum(xw)/sumw;
    header={'i' 'x_i' 'w_i' 'x_i w_i' '(x_i-M)w_i' '(x_i-M)^2 w_i'};
    dev2=((x-mx).^2).*w;
    corpus=[seq, x, w, xw, (x-mx).*w, dev2];
    sumdev2=sum(dev2);
    footer=[NaN NaN sum(w) sum(xw) 0 sumdev2];
    strtitle='Details of weighted variance calculation';
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

var=sumdev2/den;
if nargin<2
strfin=[' \it VAR$(X)= \frac{\sum_{i=1}^n (x_i -M)^2}{n-1}= \frac{'  num2str(sumdev2) '}{' num2str(den) '}=' num2str(var) '$'];
else
strfin=[' \it VAR$(X)= \frac{\sum_{i=1}^n (x_i -M)^2 w_i}{\sum_{i=1}^n w_i}= \frac{'  num2str(sumdev2) '}{' num2str(den) '}=' num2str(var) '$'];  
end

fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strfin,'Interpreter','latex','FontSize',fs1);

end
