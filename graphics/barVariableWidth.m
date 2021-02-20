function barVariableWidth(heights, classes, varargin)
%barVariableWidth produces a bar plot with different widths and colors for each bar
%
%<a href="matlab: docsearchFS('barVariableWidth')">Link to the help function</a>
%
% Required input arguments:
%
%  heights :  y-coordinates. Vector. Vector
%              Vector of length k containing numeric values describing the
%              bars which make up the plot.
%
%   classes :  classes of the frequency distribution. Vector.
%              Vector of length k+1 containing numeric values describing the
%              widths of the bars which make up the plot. classes(1)
%              contains the starting point of the first bar on the x axis, and
%              classes(end) contains the end point of the last bar on the x axis.
%              For example if classes=[0.5 0.6 0.9 1 1.2]; the first baar
%              has x coordinates [0.5 0.6], the second bar has x
%              coordinates [0.6 0.9], the third has x coordinates [0.9 1]
%              and the fourth has x coordinates [1 1.2]
%
% Optional input arguments:
%
%         Color :   Rectangle colors. scalar | vector | matrix | RGB triplet | 'r' | 'g' | 'b' | .
%                   Rectangle colors, specified as a scalar, vector,
%                   matrix, or a color name.
%                   Example - 'color',1:5
%                   Data Types - scalar | vector | matrix | RGB triplet | 'r' | 'g' | 'b' | .
%
%    LineWidth :   Line Width of the vertices. Scalar.
%                   Scalar containing the width of the lines of the
%                   rectangles.
%                   Example - 'LineWidth',2
%                   Data Types - double
%
%
%    FaceAlpha :   Face transparency.
%                   1 (default) | scalar in range [0,1] | 'flat' | 'interp'.
%                   A value of 1 is fully opaque and 0 is completely transparent.
%                   For additional details about this option see option
%                   FaceAlpha inside patch.
%                   Example - 'FaceAlpha',0.8
%                   Data Types - double
%
%    EdgeColor :   Edge colors.
%                   [0 0 0] (default) | 'none' | 'flat' | 'interp' | RGB triplet | hexadecimal color code | 'r' | 'g' | 'b' | .
%                   Colors of the edges of the rectangles.
%                   For additional details about this option see option
%                   FaceAlpha inside patch.
%                   Example - 'EdgeColor',[0 0.5 1]
%                   Data Types -  | 'none' | 'flat' | 'interp' | RGB triplet | hexadecimal color code | 'r' | 'g' | 'b' |
%
%    LineStyle :   Line style.
%                   '-' (default) | '--' | ':' | '-.' | 'none''
%                   Line Style of the edges of the rectangles.
%                   Example - 'LineStyle','--']
%                   Data Types -  '-' (default) | '--' | ':' | '-.' | 'none''
%
%  Output:
%
%
% See also: bar
%
% References:
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('barVariableWidth')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % barVariableWidth with all default options.
    close all
    % The following table shows the frequency distribution of the firms in
    % correspondence of the different classes of number of employees
    labels={'10-20' '20-30' '30-50' '50-100' '100-170'};
    freqs=[132; 180; 60; 48; 20];
    T=array2table(freqs,'RowNames',labels,'VariableNames',{'Frequency distribution'});
    disp(T);
    % Show the plot of the frequency densities (frequency/class width)
    % widths= vector which contains the width of the classes
    widths=[10; 10; 20; 50; 70];
    classes=[10; 20; 30; 50; 100; 170];
    dens=freqs./widths;
    % Show the plot
    barVariableWidth(dens,classes)
    xlabel('Employees')
    ylabel('Frequency density')
%}
%
%{
    %%Example in stackoverflow.
    % https://stackoverflow.com/questions/18419339/how-to-plot-bar-with-different-height-and-differenth-width-in-matlab
    close all
    x = [0.5 0.6 0.9 1 1.2 1.8];
    dy = [1 3 2 .5 .1];
    barVariableWidth(dy,x)
    xlabel('Time')
    ylabel('Prob density')
%}

%{
    %%Example in stackoverflow (different colors).
    % https://stackoverflow.com/questions/18419339/how-to-plot-bar-with-different-height-and-differenth-width-in-matlab
    close all
    x = [0.5 0.6 0.9 1 1.2 1.8];
    dy = [1 3 2 .5 .1];
    barVariableWidth(dy,x,'Color',1:5)
    xlabel('Time')
    ylabel('Prob density')
%}
%
%
%{
    %%Example in stackoverflow (option  FaceAlpha).
    % https://stackoverflow.com/questions/18419339/how-to-plot-bar-with-different-height-and-differenth-width-in-matlab
    close all
    x = [0.5 0.6 0.9 1 1.2 1.8];
    dy = [1 3 2 .5 .1];
    barVariableWidth(dy,x,'FaceAlpha',0.1)
    xlabel('Time')
    ylabel('Prob density')
    title('Option ''FaceAlpha''')
%}
%
%{
    %%Example in stackoverflow (option  LineWidth).
    % https://stackoverflow.com/questions/18419339/how-to-plot-bar-with-different-height-and-differenth-width-in-matlab
    close all
    x = [0.5 0.6 0.9 1 1.2 1.8];
    dy = [1 3 2 .5 .1];
    barVariableWidth(dy,x,'LineWidth',3)
    xlabel('Time')
    ylabel('Prob density')
    title('Option ''LineWidth''')
%}
%{
    %%Example in stackoverflow (option  EdgeColor).
    % https://stackoverflow.com/questions/18419339/how-to-plot-bar-with-different-height-and-differenth-width-in-matlab
    close all
    x = [0.5 0.6 0.9 1 1.2 1.8];
    dy = [1 3 2 .5 .1];
    barVariableWidth(dy,x,'EdgeColor','r')
    xlabel('Time')
    ylabel('Prob density')
    title('Option ''EdgeColor''')
%}

%{
    %%Example in stackoverflow (option  LineStyle).
    % https://stackoverflow.com/questions/18419339/how-to-plot-bar-with-different-height-and-differenth-width-in-matlab
    close all
    x = [0.5 0.6 0.9 1 1.2 1.8];
    dy = [1 3 2 .5 .1];
    barVariableWidth(dy,x,'LineStyle','--')
    xlabel('Time')
    ylabel('Prob density')
    title('Option ''LineStyle''')
%}

%% Beginning of code


Color='c';
LineWidth=1;
FaceAlpha=1;
EdgeColor=zeros(1,3);
LineStyle='-';

if nargin>2
    options=struct('Color',Color,'LineWidth',LineWidth,'FaceAlpha',FaceAlpha,...
        'EdgeColor',EdgeColor,'LineStyle',LineStyle);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:barVariableWidth:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    Color=options.Color;
    LineWidth=options.LineWidth;
    FaceAlpha=options.FaceAlpha;
    EdgeColor=options.EdgeColor;
    LineStyle=options.LineStyle;
end

if length(classes)-1~=length(heights)
    error('FSDA:barVariableWidth:WrongInputOpt','length(classes)-1 must be equal to length(heights).');
end

classes=classes(:)';
heights=heights(:)';

x1=classes(1:end-1);
x2=classes(2:end);
n=length(heights);

X=[repmat(x1,2,1); repmat(x2,2,1)];
zer=zeros(1,n);
Freq=[zer; repmat(heights,2,1); zer];
patch(X,Freq,Color,'LineWidth',LineWidth,'FaceAlpha',FaceAlpha,...
    'EdgeColor',EdgeColor,'LineStyle',LineStyle)

end
%FScategory:VIS-Reg