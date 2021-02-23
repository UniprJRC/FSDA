function out = GUImad(x,flag,w)
%GUImad shows the necessary calculations to obtain MAD, S_M or S_Me in a GUI.
%
%<a href="matlab: docsearchFS('GUImad')">Link to the help function</a>
%
% This routine shows all the intermediate necessary steps to compute the
% three following variability indexes:
%  \[
%   MAD =Me(|x_{i} - Me|n_i).
%  \]
%  \[
%   S_{M} = \frac{\sum_{i=1}^{r}|x_{i}-M| n_i}{n}
%  \]
%  \[
%   S_{Me} = \frac{\sum_{i=1}^r |x_i-Me|n_i}{n}.
%  \]
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
%    flag  : median or mean absolute deviation from median or mean absolute deviation from mean. Scalar.
%         If flag=1 (default),  MAD is based on medians, i.e.
%         median(abs(x-median(x)).
%         elseif flag=0, $S_M$ is computed (mean absolute deviation) i.e.
%         mean(abs(x-mean(x)).
%         elseif flag=2, $S_{Me}$ is computed (mean absolute deviation from median) i.e.
%         mean(abs(x-median(x)).
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
% See also: GUIstd, GUIvar, GUItrimmean
%
% References:
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('GUImad')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%

%{
    % Example of calculation of MAD.
    % MAD = median absolute deviation from median.
    x=[98 105 85 110 102];
    y=GUImad(x)
%}

%{
    % Example of calculation of SMe.
    % SMe= mean absolute deviation from the median.
    x=[98 105 85 110 102];
    y=GUImad(x,2);
%}

%{
    % Example of calculation of SM.
    % SM= mean absolute deviation from mean.
    x=[98 105 85 110 102];
    y=GUImad(x,0);
%}

%{
    %% MAD in a frequency distribution.
    % MAD = median absolute deviation from median.
    % Frequency distribution of the number of children in a sample of 200
    % families.
    X=[0	112
    1	156
    2	111
    3	16
    4	4
    7	1];
    x=X(:,1);
    freq=X(:,2);
    flag=1;
    GUImad(x,flag,freq);
%}

%{
    %% SM in a frequency distribution.
    % SM= mean absolute deviation from mean.
    % Frequency distribution of the number of children in a sample of 200
    % families.
    X=[0	112
    1	156
    2	111
    3	16
    4	4
    7	1];
    x=X(:,1);
    freq=X(:,2);
    flag=0;
    GUImad(x,flag,freq);
%}

%{
    %% SMe in a frequency distribution.
    % SMe= median absolute deviation from mean.
    % Frequency distribution of the number of children in a sample of 200
    % families.
    X=[0	112
    1	156
    2	111
    3	16
    4	4
    7	1];
    x=X(:,1);
    freq=X(:,2);
    flag=2;
    GUImad(x,flag,freq);
%}

%% Beginning of code

x=x(:);
lenx=length(x);

n=length(x);
seq=(1:n)';

if nargin < 2 || isempty(flag)
    flag = 1;
end

if nargin<3 % no weights
    xord=sort(x);
    Me=median(x);
    xminusMe=xord-Me;
    absxminusMe=abs(xminusMe);
    
    M=mean(x);
    xminusM=x-M;
    absxminusM=abs(xminusM);
    
    if  flag ==1 % median absolute deviation from median (MAD)
        xminusMeord=sort(absxminusMe);
        strtitle='Details of median absolute deviation calculation';
        header={'i' 'x_i' 'x_{(i)}' 'x_i-Me' '|x_i-Me|_{(i)}'};
        corpus=[seq, x, xord xminusMe xminusMeord];
        footer=[NaN sum(x) sum(xord) sum(xminusMe) sum(xminusMeord)];
        
    elseif flag==0 % mean absolute deviation from mean
        strtitle='Details of mean absolute deviation from mean ($S_M$) calculation';
        header={'i' 'x_i' 'x_i-M'  '|x_i-M|' };
        corpus=[seq, x, xminusM absxminusM];
        num=sum(absxminusM);
        footer=[NaN sum(x) 0 num];
        
        
    elseif flag==2 % mean absolute deviation from median
        strtitle='Details of mean absolute deviation from Median ($S_{Me}$)  calculation';
        header={'i' 'x_i' 'x_{(i)}' 'x_i-Me'  '|x_i-Me|' };
        corpus=[seq, x, xord xminusMe absxminusMe];
        num=sum(absxminusMe);
        footer=[NaN sum(x) sum(xord) sum(xminusMe) num];
    else
        error(message('FSDA:GUImad:BadFlagReduction'));
    end
    dim = [.15 .80 0.1 0.1];
    
else % MAD with weights
    n=sum(w);
    freq=w;
    [xord,xordind]=sort(x);
    freqxord=freq(xordind);
    cfreqord=cumsum(freqxord);
    if mod(n,2)==0
        poslow=find(cfreqord>=(n/2),1);
        posup=find(cfreqord>=(n/2)+1,1);
        xlow=xord(poslow);
        xup=xord(posup);
        Me=(xlow+xup)/2;
    else
        posMe=find(cfreqord>=(n+1)/2,1);
        Me=xord(posMe);
    end
    xminusMe=xord-Me;
    absxminusMe=abs(xminusMe);
    % end of median computation
    
    if  flag ==1 % median absolute deviation from median (MAD)
        [xminusMeord,indxminusMeord]=sort(absxminusMe);
        freqxminusMeord=freq(indxminusMeord);
        cfreqorxminusMe=cumsum(freqxminusMeord);
        
        if mod(n,2)==0
            poslowf=find(cfreqorxminusMe>=(n/2),1);
            posupf=find(cfreqorxminusMe>=(n/2)+1,1);
            xlowf=xminusMeord(poslowf);
            xupf=xminusMeord(posupf);
            MADf=(xlowf+xupf)/2;
        else
            posMe=find(cfreqorxminusMe>=(n+1)/2,1);
            MADf=xminusMeord(posMe);
        end
        
        strtitle='Details of median absolute deviation (MAD) calculation';
        header={'i' 'x_i' 'n_i' 'n_i''' 'x_i-Me' '|x_i-Me|_{(i)}' 'n_{i (|x_i-Me|)_{(i)}}' 'n_{i (|x_i-Me|)_{(i)}}''' };
        corpus=[seq, xord,   freqxord   cfreqord  xminusMe xminusMeord freqxminusMeord cfreqorxminusMe];
        footer=[NaN  sum(xord) sum(freqxord) NaN sum(xminusMe) sum(xminusMeord) sum(freqxminusMeord)  NaN];
        
    elseif flag==0 % mean absolute deviation from mean (S_M)
        strtitle='Details of mean absolute deviation from mean ($S_{M}$) calculation';
        xfreq=x.*freq;
        numM=sum(xfreq);
        M=numM/n;
        xminusM=x-M;
        absxminusM=abs(xminusM);
        absxminusMfreq=absxminusM.*freq;
        num=sum(absxminusMfreq);
        header={'i' 'x_i'  'n_i' 'x_in_i' 'x_i-M'  '|x_i-M|' '|x_i-M|n_i'};
        corpus=[seq,  x        freq  xfreq xminusM absxminusM      absxminusMfreq];
        footer=[NaN sum(x) sum(freq) numM   0      sum(absxminusM) num];
        
    elseif flag==2 % mean absolute deviation from median
        strtitle='Details of mean absolute deviation from Median ($S_{Me}$) calculation';
        header={'i' 'x_i' 'n_i' 'n_i''' 'x_i-Me'  '|x_i-Me|' '|x_i-Me|n_i'};
        absxminusMefreq=absxminusMe.*freq;
        num=sum(absxminusMefreq);
        corpus=[seq, xord freqxord cfreqord xminusMe absxminusMe absxminusMefreq];
        footer=[NaN sum(x) sum(freq) NaN sum(xminusMe) sum(absxminusMe) num];
    else
        error(message('FSDA:GUImad:BadFlagReduction'));
    end
    dim = [.02 .80 0.1 0.1];
end

str=strForSchool(header, corpus, footer);

out=array2table([corpus;footer],'VariableNames',header);

fs=14;
figure('Position',[100 100 1000 600],'Units','normalized');
annotation('textbox',dim,'FitBoxToText','on','String',str,'Interpreter','latex','FontSize',fs);

dim = [.7 .7 0.1 0.1];
if nargin<3 % no weights
    if flag==1 || flag ==2
        if mod(lenx,2)==0 % n is even
            poslow=lenx/2;
            posup=lenx/2+1;
            xlow=xord(poslow);
            xup=xord(posup);
            strmedian=['\boldmath{$Me$}= $\frac{x_{(' num2str(poslow) ')}+x_{(' num2str(posup) ')}}{' num2str(2) '}= \frac{' num2str(xlow) '+' num2str(xup) '}{' num2str(2) '}=' num2str(Me) '$'];
        else
            posMe=(lenx+1)/2;
            strmedian=['\boldmath{$Me$}= $x_{(' num2str(posMe) ')}=' num2str(Me) '$'];
        end
    elseif flag==0 % S_M
        sumx=sum(x);
        strmedian=['\boldmath{$M$}=$\frac{\sum_{i=1}^n x_i}{n}=\frac{' num2str(sumx) '}{' num2str(n) '}='  num2str(M) '$'];
    end
    
    if  flag ==1 % unweighted MAD
        if mod(lenx,2)==0 % n is even
            xlowf=xminusMeord(poslow);
            xupf=xminusMeord(posup);
            MADf=0.5*(xlowf+xupf);
            strfin=[' \it MAD$(X)= Me|x_i -Me| = \frac{ |x_i-Me|_{(n/2)}+|x_i-Me|_{(n/2+1)}}{2}=\frac{'  num2str(xlowf) '+' num2str(xupf) '}{' num2str(2) '}=' num2str(MADf) '$'];
        else % n is odd
            MADf=xminusMeord(posMe);
            strfin=[' \it MAD$(X)= Me|x_i -Me| =  |x_i-Me|_{((n+1)/2)}=|x_i-Me|_{(' num2str(posMe) ' )}='  num2str(MADf) '$'];
        end
    elseif flag==2 % S_Me unweighted mean absolute deviation from median
        MADf=num/n;
        strfin=[' \it S$_{Me}(X)= \frac{\sum_{i=1}^n |x_i -Me|}{n}= \frac{'  num2str(num) '}{' num2str(n) '}=' num2str(MADf) '$'];
        
    elseif flag==0 % S_M mean absolute deviation
        MADf=num/n;
        strfin=[' \it S$_{M}(X)= \frac{\sum_{i=1}^n |x_i -M|}{n}= \frac{'  num2str(num) '}{' num2str(n) '}=' num2str(MADf) '$'];
    else
        error(message('FSDA:GUImad:BadFlagReduction'));
    end
    
else % frequency distribution
    if  flag ==1 % MAD (median absolute deviation from median)
        if mod(n,2)==0 % n is even
            strmedian=['\boldmath{$Me$}= $\frac{x_{(' num2str(n/2) ')}+x_{(' num2str(n/2+1) ')}}{' num2str(2) '}= \frac{' num2str(xlow) '+' num2str(xup) '}{' num2str(2) '}=' num2str(Me) '$'];
            strfin=[' \it MAD$(X)= Me|x_i -Me| = \frac{ |x_i-Me|_{(n/2)}+|x_i-Me|_{(n/2+1)}}{2}=\frac{'  num2str(xlowf) '+' num2str(xupf) '}{' num2str(2) '}=' num2str(MADf) '$'];
        else
            strmedian=['\boldmath{$Me$}= $x_{(' num2str(posMe) ')}=' num2str(Me) '$'];
            strfin=[' \it MAD$(X)= Me|x_i -Me| =  |x_i-Me|_{((n+1)/2)}=|x_i-Me|_{(' num2str(posMe) ' )}='  num2str(MADf) '$'];
        end
        
    elseif flag==2  % S_Me Mean absolute deviation from median
        if mod(n,2)==0 % n is even
            strmedian=['\boldmath{$Me$}= $\frac{x_{(' num2str(n/2) ')}+x_{(' num2str(n/2+1) ')}}{' num2str(2) '}= \frac{' num2str(xlow) '+' num2str(xup) '}{' num2str(2) '}=' num2str(Me) '$'];
        else
            strmedian=['\boldmath{$Me$}= $x_{(' num2str(posMe) ')}=' num2str(Me) '$'];
        end
        MADf=num/n;
        strfin=[' \it S$_{Me}(X)= \frac{\sum_{i=1}^n |x_i -Me|n_i}{n}= \frac{'  num2str(num) '}{' num2str(n) '}=' num2str(MADf) '$'];
        
    elseif flag==0  % Mean absolute deviation from mean (S_M)
        strmedian=['\boldmath{$M$}= $\frac{\sum_{i=1}^r x_i n_i}{n}=  \frac{'  num2str(numM) '}{' num2str(n) '}$'];
        MADf=num/n;
        strfin=[' \it S$_{M}(X)= \frac{\sum_{i=1}^r |x_i -M|n_i}{n}= \frac{'  num2str(num) '}{' num2str(n) '}=' num2str(MADf) '$'];
        
    else
        error(message('FSDA:GUImad:BadFlagReduction'));
    end
end
annotation('textbox',dim,'FitBoxToText','on','String',strmedian,'Interpreter','latex','FontSize',fs);

dim = [.15 .9 0.1 0.1];
fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strtitle,'Interpreter','latex','FontSize',fs1);

% Final formula to compute the index
dim = [.2 .05 0.1 0.1];
annotation('textbox',dim,'FitBoxToText','on','String',strfin,'Interpreter','latex','FontSize',fs1);

end
%FScategory:GUI