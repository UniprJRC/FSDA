function out = GUIconcentration(Q, varargin)
%GUIconcentration shows the necessary calculations to obtain the GINI concentration index in a GUI.
%
%
%<a href="matlab: docsearchFS('GUIconcentration')">Link to the help function</a>
%
%  Required input arguments:
%
%     Q : vector of non negative numeric data representing quantities. Vector.
%           Vector containing strictly numerical data
%           Data Types - double
%
%
%  Optional input arguments:
%
%    freq : frequencies. Vector.
%         Vector of the same legth of Q containing the frequencies
%         associated to each value of x. If freq is not it is assumed that
%         all observations have the same (relative) frequency.
%           Example - 'freq',1:10
%           Data Types - double
%
% ExactFormula : use of exact formula. Boolean.
%               Use of exact or approximate formula.
%               The default is to use the exact formula if previously
%               option freq has not been supplied, else approximate formula
%               is used. Option ExactFormula overwrites the default.
%           Example - 'ExactFormula',true
%           Data Types - boolean
%
%     plots    : show Lorenz curve. Boolean.
%                If plots is true an additional plot which shows
%                the concentration area is displayed on the screen. The
%                default value of plots is false.
%           Example - 'plots',true
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
% See also: GUIvar, GUIstd, GUImad
%
% References:
%
% Cerioli, A., Milioli, M.A., Riani, M. (2016), "Esercizi di statistica (Quinta edizione)". [CMR]
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('GUIconcentration')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%
%{
    % Concentration index using exact formula.
    % The percentage breakdown of the investments of a mutual 
    % is the following.
    % This is exercise 1.29 from [CMR] - page 61.
    % Government bonds 55 per cent.
    % Bonds 12 per cent.
    % Italian equities 23 per cent
    % Foreign equities 6 per cent.
    % Cash 4 per cent.
    % Compute concentration ratio and show the concentration area.
    Q=[55 12 23 6 4];
    Q=sort(Q);
    out=GUIconcentration(Q,'plots',1);
%}

%{
    %% Concentrartion index using approximate formula.
    % The families in a certain area have been classified 
    % according to their annual income.
    % This is exercise 1.28 from [CMR] - page 59.
    % Classes of income and frequencies.
    % 5  - 10	10 
    % 10 - 20	25
    % 20 - 30	45
    % 30 - 50	65
    % 50 - 70	35
    % 70 - 100	20
    % Compute the concentration ratio and show the concentration area.
    % Note that given that n is large we can use the approximate formula.
    % Using as x_i the central value of the classes, we estimate the Q_i as
    % x_i * n_i. For example Q_1 = 7.5 * 10, Q_2 = 15 * 25...Q_6 = 85 * 20
    Q=[75 375 1125 2600 2100 1700];
    ni=[10 25 45 65 35 20];
    out=GUIconcentration(Q,'freq',ni,'plots',1);
%}

%% Beginning of code

freq='';
ExactFormula='';
plots=false;
if nargin > 1
    options=struct('freq',freq,'ExactFormula',ExactFormula,...
        'plots',plots);
    
    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:GUIconcentration:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    freq = options.freq;
    ExactFormula = options.ExactFormula;
    %         if min(freq)<0 || max(freq)>=1
    %             error('FSDA:percclass:WrongInputOpt','Optional argument prt must be in the interval [0 1].');
    %         end
    plots=options.plots;
end


n=length(Q);

if isempty(freq)
    freq=(1/n)*ones(n,1);
    
    if isempty(ExactFormula)
        approx=false;
    end
else
    freq=freq(:)/sum(freq);
    
    if isempty(ExactFormula)
        approx=true;
    end
end


if ~isempty(ExactFormula)
    approx=ExactFormula;
end

Q=Q(:);
seq=(1:n)';
if approx==true % approximate formula is used
    header={'i' 'Q_{i}' 'f_i' 'f_i''' 'q_i' 'q_i''' '(q_i''+q_{i-1}'')f_i'};
    
    freqcum=cumsum(freq);
    q=Q/sum(Q);
    qcum=cumsum(q);
    qcum2=[qcum(1);  qcum(2:n)+qcum(1:n-1)].*freq;
    corpus=[seq, Q freq, freqcum, q, qcum, qcum2];
    sumqcum=sum(qcum);
    sumqcum2=sum(qcum2);
    
    footer=[NaN NaN 1 NaN  1 sumqcum sumqcum2];
    strtitle='Details of Gini concentration index calculation (approx formula)';
else % exact formula is used
    header={'i' 'Q_{i}' 'f_i' 'f_i''' 'q_i' 'q_i''' '(q_i''+q_{i-1}'')f_i'};
    
    freqcum=cumsum(freq);
    q=Q/sum(Q);
    qcum=cumsum(q);
    qcum2=[qcum(1);  qcum(2:n)+qcum(1:n-1)].*freq;
    corpus=[seq, Q, freq, freqcum, q, qcum, qcum2];
    sumqcum=sum(qcum);
    sumqcum2=sum(qcum2);
    
    footer=[NaN NaN 1 NaN  1 sumqcum sumqcum2];
    strtitle='Details of Gini concentration index calculation (exact formula)';
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

dim = [.2 .9 0.1 0.1];
fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strtitle,'Interpreter','latex','FontSize',fs1);

dim = [.2 .05 0.1 0.1];

out=struct;

if approx==true
    sumf=footer(end);
    Q=1- sumf;
    strfin=[' \it R $\approx 1- \sum_{i=1}^r (q_i +q_{i-1}) f_i = 1 -'  num2str(sumf) '= ' num2str(Q) '$'];
    out.R=Q;
else
    % Write details to find final calculation
    finalR=((n+1)/(n-1))-(2/(n-1))*sumqcum;
    strfin=[' \it R $= \frac{n+1}{n-1} - \frac{2}{n-1} \sum_{i=1}^r q_i''  = \frac{' num2str(n+1) ...
        '}{'  num2str(n-1) '} -\frac{2}{'  num2str(n-1) '}' num2str(sumqcum) '='  num2str(finalR) '$'];
    out.R=finalR;
end

fs1=20;
annotation('textbox',dim,'FitBoxToText','on','String',strfin,'Interpreter','latex','FontSize',fs1);

if plots==true
    figure
    % Show graphically concentration area
    f1 = [0; freqcum];
    q1 = [0; qcum];
    fill([0;f1],[0;q1],'c');
    
    hold on
    % plot([0;f1],[0;q1],'o','Marker','o');
    plot([0,1],[0,1],'--k');     
    stem(f1,q1,'--')                
    axis tight      % ranges of abscissa and ordinate are exactly [0,1]
    axis square     % same length on hor and ver axes
    set(gca,'XTick',get(gca,'YTick'))   % ensure equal ticking
    xlabel('$f_i''$','Interpreter','latex','FontSize',fs1)
    ylabel('$q_i''$','Interpreter','latex','FontSize',fs1)
end

out.data=array2table([corpus;footer],'VariableNames',header);

end
%FScategory:GUI