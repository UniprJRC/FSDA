function out = GUIquantile(x, z, varargin)
%GUIquantile shows the necessary calculations to obtain $x_z$ quantile.
%
%
%<a href="matlab: docsearchFS('GUIquantile')">Link to the help function</a>
%
%  Required input arguments:
%
%     x : vector of numeric data. Vector.
%           Vector containing strictly numerical data.
%           Note that if x is referred to a continuous variable x(1)
%           represent the lower extreme of the first class, and the other
%           values of x contain the upper extremes of the other
%           classes.
%           Data Types - double, ordered categorical
%
%     z : requested probability. Scalar.
%           Requested probabilities for which to compute the quantile,
%           specified as a scalar.
%           Data Types - double
%
%  Optional input arguments:
%
%    freq : frequencies. Vector.
%         Vector of the same length of x containing the frequencies
%         (weights) associated to each value of x if x is discrete or
%         vector with length (x)-1 if variable X is continuous. If
%         freq is not supplied, it is assumed that all observations have the
%         same (relative) frequency (weight).
%           Example - 'freq',1:10
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
%     plots    : show quantile graphically. Boolean.
%                If plots is true an additional plot which shows
%                graphically how the quantile has been obtained using the
%                linear interpolation is displayed on the screen. The
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
% See also: GUIconcentration, GUIvar, GUIstd
%
% References:
% Cerioli, A., Milioli, M.A., Riani, M. (2016), "Esercizi di statistica (Quinta edizione)". [CMR]
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('GUIquantile')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%
%{
    % Percentile "z" in an individual series.
    % In this case we find percentile 43, all frequencies are equal to 1. 
    x=[2 10 5 9 13];
    GUIquantile(x,0.43)
%}

%{
    % Computation of the median in presence of discrete data with frequencies.
    % X = Number of components of Italian families (source ISTAT).
    % freq = Number of families (in thousands).
    x=1:6;
    freq=[7910 6833 5116 4051 1088  303];
    GUIquantile(x,0.5,'freq',freq,'DiscreteData',true)
%}


%{
    % Computation of the median in presence of discrete data with
    % frequencies and related plot. (See page 14 of [CMR])
    % X = 133 students grades on a given exam
    % freq = frequencies.
    x=18:30;
    freq=[12 10 4 8 7 19 22 8 9 13 10 8 3 ];
    % Find 83 percentile and show associated plot.
    GUIquantile(x,0.83,'freq',freq,'DiscreteData',true,'plots',true)
%}

%{
    %% Example of computation of 40 per cent percentile in a frequency
    % distribution (X is continuous).
    % The following frequency distribution shows the amount (in thousands
    % of Euros) of advertising expenditure made in a given month by a series of
    % companies. (See page 15 of [CMR])
    X=[9  1787
    15  1310
    19  972
    25  2753
    35  4227
    50  2174
    100 920
    150  138
    300  54
    Inf  9];
    x=[0;X(:,1)];
    freq=X(:,2);
    GUIquantile(x,0.42,'freq',freq,'DiscreteData',false,'plots',1)
%}

%% Beginning of code

freq='';
DiscreteData=true;
plots=false;
if nargin > 2
    options=struct('freq',freq,'DiscreteData',DiscreteData,...
        'plots',plots);

    UserOptions=varargin(1:2:length(varargin));
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:GUIquantile:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    freq = options.freq;
    %         if min(freq)<0 || max(freq)>=1
    %             error('FSDA:percclass:WrongInputOpt','Optional argument prt must be in the interval [0 1].');
    %         end
    plots=options.plots;
    DiscreteData=options.DiscreteData;
end


if isempty(freq)
    sorx=sort(x(:));
    n=length(x);
    prob=(0.5:n-0.5)'/n;
    close all
    xadd=[min(x);sorx];
    probadd=[0;prob];
    figure
    % Make sure that that figure is also visible inside .mlx files
    scatter([],[]);
    axis('off')
    % set(gcf,'Visible','on')
    plot(xadd,probadd,'r','Marker','o')

    xlabel('x and requested quantile')
    ylabel('Probability')
    % plot(probadd,xadd,'r')
    x=double(x);
    quan=quantile(x,z);
    % interp1(probadd,xadd,0.4)
    hold('on')
    plot([quan; quan; min(x)],[0; z; z],'k--')
    text(quan, 0.05,['x_{' num2str(z) '}=' num2str(quan)],'FontSize',16)
    % freqcum=cumsum((1/n)*ones(n,1))
    h=cdfplot(x);
    h.Color='b';
    pause(0.0001)
    clickableMultiLegend({'Sorted values and associated quantiles'  ...
        ['Requested quantile corresponding to p=' num2str(z)] 'F(x)'},'Location','northwest');
    out=struct;
    %out.data=array2table([corpus;footer],'VariableNames',header);
    out.quantile=quan;

else
    if DiscreteData == false
        [~,indminx]=min(x(:)) ;
        xori=x;
        x(indminx)=[];
    end

    [x,indsorx]=sort(x(:));
    k=length(x);

    freq=freq(indsorx(1:end));

    freq=freq(:);
    n=sum(freq);
    seq=(1:k)';
    header={'i' 'x_{(i)}' 'n_{i}' 'f_i' 'F(x_i)' };
    f=freq/n; % relative frequencies
    fcum=cumsum(f); % cumulative relative frequencies


    corpus=[seq, x, freq, f, fcum];

    footer=[NaN NaN sum(freq) 1 NaN];

    if DiscreteData == true
        strtitle='Details of quantile computation (X is discrete)';
        classes='';
    else % continuous data
        strtitle='Details of quantile computation (X is continuous)';
        classes=strcat(string(xori(1:end-1)),'-',string(xori(2:end)));
    end


    str=strForSchool(header, corpus, footer,classes);



    fs=14;
    dim = [.2 .80 0.1 0.1];
    figure('Position',[100 100 1100 600],'Units','normalized');
    % Make sure that that figure is also visible inside .mlx files
    scatter([],[]);
    axis('off')
    % set(gcf,'Visible','on')
    annotation('textbox',dim,'FitBoxToText','on','String',str,'Interpreter','latex','FontSize',fs);

    dim = [.2 .9 0.1 0.1];
    fs1=20;
    annotation('textbox',dim,'FitBoxToText','on','String',strtitle,'Interpreter','latex','FontSize',fs1);


    if DiscreteData == true
        dim = [.2 .05 0.1 0.1];
        zstr=num2str(z);
        if z==0.5 && mod(n,2)==0
            indexxz1=find(fcum>=z,1,'first');
            indexxz2=find(fcum>z,1,'first');
            xz=0.5*(x(indexxz1)+x(indexxz2));

            xzstr=num2str(xz);
            strfin=[' \it x $_{'  zstr '}= \frac{' num2str(x(indexxz1)) '+' num2str(x(indexxz2)) '}{2}=' num2str(xz) '$'];

        else
            indexxz=find(fcum>=z,1,'first');
            xz=x(indexxz);
            xzstr=num2str(xz);
            strfin=[' \it x $_{'  zstr '}= ' num2str(xz) '$'];

        end
    else
        dim = [.02 .05 0.1 0.1];
        indexxz=find(fcum<=z,1,'last');
        zstr=num2str(z);
        xbars=x(indexxz);
        xzstr=num2str(xbars);
        reqq=xbars+(x(indexxz+1)-xbars)/(fcum(indexxz+1)-fcum(indexxz))*(z-fcum(indexxz));
        % Write details to find final calculation
        calc=[xzstr '+ \frac{' num2str(x(indexxz+1)) '-'...
            num2str(x(indexxz)) '}{' num2str(fcum(indexxz+1)) '-' num2str(fcum(indexxz)) '}'...
            '(' zstr '-' num2str(fcum(indexxz)) ')=' num2str(reqq) ]   ;
        strfin=[' \it x $_{' zstr '}= \overline x_s+ \frac{\overline{\overline x}_s  -\overline x_s}{F(x_s)-F(x_{s-1})}[' zstr '-F(x_{s-1})]=' calc '$'];
        xz=reqq;
        xzstr=num2str(reqq);
    end

    fs1=20;
    annotation('textbox',dim,'FitBoxToText','on','String',strfin,'Interpreter','latex','FontSize',fs1);

    if plots==true
        if DiscreteData == true
            figure
            bar(x,fcum,'cyan')
        else
            figure
            barVariableWidth(fcum, xori,'Color','w')
            x=xori;
        end

        hold('on')
        xlimmin=min(x)-0.5;
        xlimmax=max(x)+0.5;
        xlim([xlimmin xlimmax])

        plot([xz; xz; xlimmin],[0; z; z],'k--')
        text(xz, 0.05,['x_{' num2str(z) '}=' xzstr],'FontSize',16)
        xlabel('$x_i$','Interpreter','latex','FontSize',16)
        ylabel('Cumulative distribution function $F(x_i)$','Interpreter','latex','FontSize',16)

        if DiscreteData==false
            plot([xbars;x(indexxz+2)],fcum(indexxz:indexxz+1))
        end
    end

    out=struct;
    out.data=array2table([corpus;footer],'VariableNames',header);
    out.quantile=xz;
end

end
%FScategory:GUI