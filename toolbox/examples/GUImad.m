function out = GUImad(x, flag, w, varargin)
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
%           Vector containing strictly numerical data. Note that if x is
%           referred to a continuous variable (that is option DiscreteData
%           is false) x(1) represents the lower extreme of the first class,
%           and the other values of x contain the upper extremes of the
%           other classes.
%   Data Types - double, ordered categorical
%
%
%  Optional input arguments:
%
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
%         Vector of the same length of x containing the frequencies
%         (weights) associated to each value of x if x is discrete or
%         vector with length (x)-1 if variable X is continuous. If
%         w is not supplied, it is assumed that all observations have the
%         same (relative) frequency (weight).
%           Example - 1:10
%           Data Types - double
%
% DiscreteData : Discrete data or continuous data. Boolean.
%               If DiscreteData is true (default), we assume that the data
%               in x come from a discrete variable.
%               If DiscreteData is false, we assume that the data
%               in x come from a continuous variable.
%           Example - 'DiscreteData',false
%           Data Types - boolean
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
% Cerioli, A., Milioli, M.A., Riani, M. (2016), "Esercizi di statistica (Quinta edizione)". [CMR]
%
% Copyright 2008-2025.
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
    y=GUImad(x);
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
    %% MAD in a frequency distribution (discrete data).
    % MAD = median absolute deviation from median.
    % Frequency distribution of the number of children in a sample of 400
    % families. (See page 29 of [CMR]).
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
    %% SM in a frequency distribution (discrete data).
    % SM= mean absolute deviation from mean.
    % Frequency distribution of the number of children in a sample of 400
    % families. (See page 29 of [CMR]).
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
    %% SMe in a frequency distribution (discrete data).
    % SMe= mean absolute deviation from the median
    % Frequency distribution of the number of children in a sample of 400 families. (See page 29 of [CMR]).
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


%{
   %% MAD in a frequency distribution (continuous data).
    % Conside the following distribution of employees of a large company
    % according to monthly gross wages:
    % retr= classes of retribution (Euros)
    retr=[1000 1200; 1200 1500; 1500 2000; 2000 2500; 2500 3500; 3500 5000];
    % freq= frequencies
    freq=[30; 130; 150; 50; 30; 20];
    % x = vector which contains the extremes of the classes
    x=[retr(1,1);retr(:,2)];
    flag=1;
    % GUImad is called with input option 'DiscreteData',false
    GUImad(x,flag,freq,'DiscreteData',false);    
%}

%{
    % SMe in a frequency distribution (continuous data).
    % SMe= median absolute deviation from mean.
    % Conside the following distribution of employees of a large company
    % according to monthly gross wages:
    % retr= classes of retribution (Euros)
    retr=[1000 1200; 1200 1500; 1500 2000; 2000 2500; 2500 3500; 3500 5000];
    % freq= frequencies
    freq=[30; 130; 150; 50; 30; 20];
    % x = vector which contains the extremes of the classes
    x=[retr(1,1);retr(:,2)];
    flag=2;
    % GUImad is called with input option 'DiscreteData',false
    GUImad(x,flag,freq,'DiscreteData',false);    
%}

%% Beginning of code

x=x(:);
lenx=length(x);


if nargin<2 || isempty(flag)
    flag = 1;
end

DiscreteData=true;

if nargin > 3
    options=struct('DiscreteData',DiscreteData);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:GUImad:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    aux.chkoptions(options,UserOptions)

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    DiscreteData=options.DiscreteData;
end

seq=(1:length(x))';

if nargin<3  || isempty(w) % no weights
    n=length(x);

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
        strtitle='Details of mean absolute deviation from Median ($S_{Me}$) calculation';
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

    if DiscreteData ==true

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

    else % Continuous data
        xori=x;
        freqxord=freq;
        cfreqord=cumsum(freqxord);
        xord=0.5*(xori(1:end-1)+xori(2:end));
        Me=interp1(cumsum(freq)/n,xori(2:end),0.5);
        seq=seq(1:end-1);
        x=xord;

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


if DiscreteData == true
    classes='';
else % continuous data
    strtitle=[strtitle ' (X is continuous)'];
    classes=strcat(string(xori(1:end-1)),'-',string(xori(2:end)));
end

str=strForSchool(header, corpus, footer,classes);


fs=14;
figure('Position',[100 100 1000 600],'Units','normalized');
% Make sure that that figure is also visible inside .mlx files
scatter([],[]);
axis('off')
set(gcf,'Visible','on')
annotation('textbox',dim,'FitBoxToText','on','String',str,'Interpreter','latex','FontSize',fs);


if nargin<3 % no weights
    if flag==1 || flag ==2
        if verLessThanFS('25.1')
            indexeq='\boldmath{$Me$}= $';
        else
            indexeq= '$Me=';
        end

        if mod(lenx,2)==0 % n is even
            poslow=lenx/2;
            posup=lenx/2+1;
            xlow=xord(poslow);
            xup=xord(posup);
            strmedian=[indexeq ' \frac{x_{(' num2str(poslow) ')}+x_{(' num2str(posup) ')}}{' num2str(2) '}= \frac{' num2str(xlow) '+' num2str(xup) '}{' num2str(2) '}=' num2str(Me) '$'];

        else
            posMe=(lenx+1)/2;

            strmedian=[indexeq 'x_{(' num2str(posMe) ')}=' num2str(Me) '$'];
        end
    elseif flag==0 % S_M
        sumx=sum(x);
        if verLessThanFS('25.1')
            indexeq='\boldmath{$M$}= $';
        else
            indexeq= '$M=';
        end

        strmedian=[ indexeq '\frac{\sum_{i=1}^n x_i}{n}=\frac{' num2str(sumx) '}{' num2str(n) '}='  num2str(M) '$'];
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
        if DiscreteData == true
            if mod(n,2)==0 % n is even
                strmedian=['\boldmath{$Me$}= $\frac{x_{(' num2str(n/2) ')}+x_{(' num2str(n/2+1) ')}}{' num2str(2) '}= \frac{' num2str(xlow) '+' num2str(xup) '}{' num2str(2) '}=' num2str(Me) '$'];
                strfin=[' \it MAD$(X)= Me|x_i -Me| = \frac{ |x_i-Me|_{(n/2)}+|x_i-Me|_{(n/2+1)}}{2}=\frac{'  num2str(xlowf) '+' num2str(xupf) '}{' num2str(2) '}=' num2str(MADf) '$'];
            else
                strmedian=['\boldmath{$Me$}= $x_{(' num2str(posMe) ')}=' num2str(Me) '$'];
                strfin=[' \it MAD$(X)= Me|x_i -Me| =  |x_i-Me|_{((n+1)/2)}=|x_i-Me|_{(' num2str(posMe) ' )}='  num2str(MADf) '$'];
            end
        else
            posMe=find(xminusMeord<=MADf,1,'last');
            %      strmedian=['\boldmath{$Me$}= $x_{(' num2str(posMe) ')}=' num2str(Me) '$'];
            strmedian=['\boldmath{$Me$}=' num2str(Me) ];
            strfin=[' \it MAD$(X)= Me|x_i -Me| =  |x_i-Me|_{(' num2str(posMe) ' )}='  num2str(MADf) '$'];

        end

    elseif flag==2  % S_Me Mean absolute deviation from median
        if DiscreteData == true
            if mod(n,2)==0 % n is even
                strmedian=['\boldmath{$Me$}= $\frac{x_{(' num2str(n/2) ')}+x_{(' num2str(n/2+1) ')}}{' num2str(2) '}= \frac{' num2str(xlow) '+' num2str(xup) '}{' num2str(2) '}=' num2str(Me) '$'];
            else
                strmedian=['\boldmath{$Me$}= $x_{(' num2str(posMe) ')}=' num2str(Me) '$'];
            end
            MADf=num/n;
            strfin=[' \it S$_{Me}(X)= \frac{\sum_{i=1}^n |x_i -Me|n_i}{n}= \frac{'  num2str(num) '}{' num2str(n) '}=' num2str(MADf) '$'];
        else
            strmedian=['\boldmath{$Me$}=' num2str(Me) ];
            MADf=num/n;
            strfin=[' \it S$_{Me}(X)= \frac{\sum_{i=1}^n |x_i -Me|n_i}{n}= \frac{'  num2str(num) '}{' num2str(n) '}=' num2str(MADf) '$'];

        end

    elseif flag==0  % Mean absolute deviation from mean (S_M)

        strmedian=['\boldmath{$M$}= $\frac{\sum_{i=1}^r x_i n_i}{n}=  \frac{'  num2str(numM) '}{' num2str(n) '}$=' num2str(M) ];
        MADf=num/n;
        strfin=[' \it S$_{M}(X)= \frac{\sum_{i=1}^r |x_i -M|n_i}{n}= \frac{'  num2str(num) '}{' num2str(n) '}=' num2str(MADf) '$'];

    else
        error(message('FSDA:GUImad:BadFlagReduction'));
    end
end

if DiscreteData==true
    dimstrmedian = [.7 .7 0.1 0.1];
else
    dimstrmedian = [.8 .7 0.1 0.1];
end
annotation('textbox',dimstrmedian,'FitBoxToText','on','String',strmedian,'Interpreter','latex','FontSize',fs);

dim = [.05 .9 0.1 0.1];
fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strtitle,'Interpreter','latex','FontSize',fs1);

% Final formula to compute the index
dim = [.2 .05 0.1 0.1];
annotation('textbox',dim,'FitBoxToText','on','String',strfin,'Interpreter','latex','FontSize',fs1);

out=struct;
out.data=array2table([corpus;footer],'VariableNames',header);
out.mad=MADf;


end
%FScategory:GUI