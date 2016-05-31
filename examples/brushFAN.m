function brushFAN(eventdata)
% brushFAN displays a GUI which enables brushing in the fanplot
%
%<a href="matlab: docsearchFS('brushFAN')">Link to the help page for this function</a>
%
%brushFAN displays a GUI where it is possible to brush steps from the fan plot
% and to see the corresponding units highlighted in other plots
%
% Required input arguments:
%
% Optional input arguments:
%
% eventdata  : scalar integer (from 1 to 3). Automatic code execution
%              without user interaction. This option enables to perform in
%              an automatic way the code associated with a particular
%              radiobutton in the GUI 
%              Example - 2 (the example associated
%              with the second radiobutton will be automatically executed)
%              Data Types - integer
%
%  Output: 
%
%
% See also: brushRES, brushROB
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('brushFAN')">Link to the help page for this function</a>
% Last modified mar 17 mag 2016 12:19:52

% Examples:

%{
    % Interactive_example
    % Display a GUI where it is possible to brush steps from the fan plot
    % and to see the corresponding units highlighted in other plots
    brushFAN
%}
%
%{
    %% Run examples associated with radiobuttons 1 to 4
    for j=1:4
        brushFAN(j);
    end
%}
%
%{
    %% Run the example associated with radiobutton 2
    brushFAN(2);
%}
%
%{
    %% Run the example associated with radiobutton 3
    brushFAN(3);
%}
%

%% Beginning of code

if nargin < 1
    % other open demos create problems. Delete them before starting this new one.
    delete(((findobj('type','figure','Tag','demo'))));
    
    figure; % Create the button group.
    h = uibuttongroup('visible','off','Position',[0 0 .3 1]);
    
    % Set name of the GUI
    set(gcf,'Name', 'Example of robust dynamic brushing starting from the FAN plot', 'NumberTitle', 'off');
    
    wid=0.95;
    sx=0.05;
    hei=0.04;
    fontsizebig = 10;
    fontsizesmall = 9;
    fontname = 'Lucida';
    backgroundcolor = 'w';
    
    examp=which('examples_regression.m');
    examp1=strrep(examp,'\','\\');
    
    % Include on the left the radio buttons
    uicontrol('Style','Radio','String','Wool dataset','Units','normalized',...
        'pos',[sx 0.9 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
    uicontrol('Style','Radio','String','Stack loss dataset','Units','normalized',...
        'pos',[sx 0.8 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
    uicontrol('Style','Radio','String','Loyalty cards','Units','normalized',...
        'pos',[sx 0.7 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
    
    
    str = sprintf(['\n'...
        'To explore how dynamic brushing works, select one of the radio buttons on the left. \n\n' ...
        'The fan plot will appear automatically. \n\n'...
        'The points associated with the trajectories you select in the fan plot will also be automatically highlighted in the plots of monitoring residuals and in the scatter plot matrix.\n'...
        '\n\n'...
        'To produce the FAN plot in order to decide which is the most appropriate value of BOX-COX lambda, type:\n\n\n' ...
        '        [out]=FSRfan(y,X) \n' ...
        '        fanplot(out) \n\n\n' ...
        'See the help files of the above functions for more information \n\n'...
        'or see file']);
    
    annotation('textbox',[0.3 0 0.7 1],'String',str,'Interpreter','none','FontSize',fontsizebig,'FontName',fontname,'BackgroundColor',backgroundcolor);
    
    import com.mathworks.mlwidgets.html.HTMLRenderer;
    
    % create component
    r = HTMLRenderer;
    % set the text to display
    r.setHtmlText(['<html> <a href="matlab:opentoline(''',examp1,''',538)">example_regression.m</a></html>']);
    % make sure the component is opaque
    r.setOpaque(true);
    % add the component
    javacomponent(r, [290 85 400 35], gcf);
    %[left, bottom, width, height]
    
    % Include at the bottom of the GUI two toggle buttons named close and
    % closeall
    uicontrol('Style','Pushbutton','Units','normalized', ...
        'Position',[0.01 0.07 0.17 0.05], ...
        'Callback','close', 'String','Close');
    uicontrol('Style','Pushbutton','Units','normalized', ...
        'Position',[0.01 0.01 0.28 0.05], ...
        'Callback','cabc', 'String','Close all plots but current GUI');
    
    
    % Initialize some button group properties.
    set(h,'SelectionChangeFcn',@brushFANex);
    set(h,'SelectedObject',[]);  % No selection
    set(h,'Visible','on'); % Make the object visible
    
    % Add the logo to the GUI
    imdata = imread('logo.png','BackgroundColor',(240/255)*[1 1 1]);
    hA=axes('Position',[0.05 0.15 0.2 0.2]);
    image(imdata,'Parent',hA);
    set(hA,'DataAspectRatio',[1 1 1]);
    axis('off')
else
    brushFANex([],eventdata)
end

    function brushFANex(~,eventdata)
        
        if eventdata==1 || (isa(eventdata,'matlab.ui.eventdata.SelectionChangedData') && strcmp(get(eventdata.NewValue,'String'),'Wool dataset'));
            wool=load('wool.txt');
            y=wool(:,4);
            X=wool(:,1:3);
            [out]=FSRfan(y,X);
            namey='Number of cycles to failure';
            nameX={'Length of test specimen', 'Amplitude of loading cycle', 'Load '};
            if isa(eventdata,'matlab.ui.eventdata.SelectionChangedData')
                
                fanplot(out,'lwd',1.5,'FontSize',11,'SizeAxesNum',11,'nameX',nameX,'namey',namey,'databrush',{'selectionmode' 'Brush'...
                    'persist' '' 'multivarfit' '2' 'FlagSize' '5' 'Label' 'on' 'RemoveLabels' 'off'})
            else
                fanplot(out,'lwd',1.5,'FontSize',11,'SizeAxesNum',11,'nameX',nameX,'namey',namey,'databrush','')
            end
        elseif eventdata==2 || (isa(eventdata,'matlab.ui.eventdata.SelectionChangedData') && strcmp(get(eventdata.NewValue,'String'),'Stack loss dataset'))
            
            
            stack_loss_data=load('stack_loss.txt');
            y=stack_loss_data(:,4);
            X=stack_loss_data(:,1:3);
            namey='Number of cycles to failure';
            nameX={'Length of test specimen', 'Amplitude of loading cycle', 'Load '};
            [out]=FSRfan(y,X,'init',8);
            if isa(eventdata,'matlab.ui.eventdata.SelectionChangedData')
                
                fanplot(out,'ylimy',[-5 7],'lwd',1.5,'FontSize',11,'SizeAxesNum',11,'nameX',nameX,'namey',namey,'databrush',{'selectionmode' 'Brush'...
                    'persist' '' 'multivarfit' '2' 'FlagSize' '5' 'Label' 'on' 'RemoveLabels' 'off'})
            else
                fanplot(out,'ylimy',[-5 7],'lwd',1.5,'FontSize',11,'SizeAxesNum',11,'nameX',nameX,'namey',namey,'databrush','')
            end
            
        elseif eventdata==3 || (isa(eventdata,'matlab.ui.eventdata.SelectionChangedData') && strcmp(get(eventdata.NewValue,'String'),'Loyalty cards'))
            
            
            loyalty=load('loyalty.txt');
            y=loyalty(:,4);
            X=loyalty(:,1:3);
            % Compute fan plot to find best value of transformation parameter
            [out]=FSRfan(y,X,'la',[0 1/3 0.4 0.5]);
            % Dynamic Brushing starting from the fan plot
            %Example of the use of FlagSize, namey, namex, lwd,FontSize, SizeAxesNum.
            namey='Sales';
            nameX={'Number of visits', 'Age', 'Number of persons in the family'};
            %FlagSize controls how large must be the highlighted points. It is a
            %parameter of selectdataFS.
            if isa(eventdata,'matlab.ui.eventdata.SelectionChangedData')
                
                fanplot(out,'ylimy',[-10 20],'xlimx',[10 520],'lwd',1.5,'FontSize',11,'SizeAxesNum',11,'nameX',nameX,'namey',namey,'databrush',{'selectionmode' 'Brush'...
                    'multivarfit' '2' 'FlagSize' '5'})
            else
                fanplot(out,'ylimy',[-10 20],'xlimx',[10 520],'lwd',1.5,'FontSize',11,'SizeAxesNum',11,'nameX',nameX,'namey',namey,'databrush','')
            end
            
        else
        end
        
        % Now position the plots in particular areas of the screen
        scrsz = get(0,'ScreenSize');
        % Check which plots are open
        % ScreenSize is a four-element vector: [left, bottom, width, height]:
        width=scrsz(3)/3;
        
        % Check if figure containing mdr is present and get its handle
        plmdr=((findobj('type','figure','Tag','pl_mdr')));
        
        % Check if figure containing residulas is present and get its handle
        plresfwd=((findobj('type','figure','Tag','pl_resfwd')));
        
        % Check if yX is present and get its handle
        plyX=((findobj('type','figure','Tag','pl_yX')));
        
        % Check if fan plot is present and get its handle
        plfan=((findobj('type','figure','Tag','pl_fan')));
        
        % Set the position of the plots
        set(plfan,'Position',[10 scrsz(4)/10  width scrsz(4)/3])
        set(plresfwd,'Position',[(width+10) scrsz(4)/10  width scrsz(4)/3])
        set(plyX,'Position',[(2*width+10) scrsz(4)/10  width scrsz(4)/3])
        set(plmdr,'Position',[(3*width+10) scrsz(4)/10  width scrsz(4)/3])
        
        
        
    end

a=version;
if str2double(a(1))>=8
    stri='Detailed information about the datasets used in this GUI can be found in the <a href="matlab: doc -classic">section USER GUIDE of the FSDA html help system </a>';
else
    stri='Detailed information about the datasets used in this GUI can be found <a href="matlab: docsearchFS(''datasets_reg'')">here</a>';
end
disp(stri)

end
%FScategory:GUI