function brushROB(eventdata)
% brushRES displays a GUI which enables brushing in resindexplot
%
%<a href="matlab: docsearchFS('brushROB')">Link to the help page for this function</a>
%
%brushROB displays a GUI to brush units from the index plot of residuals
%and see the corresponding units highlighted in the yXplot
%
% Required input arguments:
%
% Optional input arguments:
%
% eventdata  : scalar integer (from 1 to 6). Automatic code execution
%              without user interaction. This option enables to perform in
%              an automatic way the code associated with a particular
%              radiobutton in the GUI 
%              Example - 2 (the example associated
%              with the second radiobutton will be automatically executed)
%              Data Types - integer
%
% See also: brushFAN
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('brushROB')">Link to the help page for this function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    % Interactive_example
    % Display a GUI to brush units from the index plot of residuals
    % and see the corresponding units highlighted in the yXplot
    brushROB
%}
%
%{
    %% Run examples associated with radiobuttons 1 to 6
    for j=1:6
        brushROB(j);
    end
%}
%

%% Beginning of code

if nargin < 1
    % other open demos create problems. Delete them before starting this new one.
    delete(((findobj('type','figure','Tag','demo'))));
    
    figure; % Create the GUI figure and set name of the GUI
    set(gcf,'Name', 'Example of robust dynamic brushing starting from the index plot of residuals', 'NumberTitle', 'off','Tag','demo');
    
    % Create the button group
    h = uibuttongroup('visible','off','Position',[0 0 .3 1]);
    
    % Size of the radio buttons
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
    
    uicontrol('Style','Radio','String','Forbes dataset','Units','normalized',...
        'pos',[sx 0.9 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
    annotation('textbox',[sx 0.87 0.1 hei],'String','LMS','EdgeColor','none','FontSize',fontsizesmall,'FontName',fontname);
    %
    uicontrol('Style','Radio','String','Multiple regression','Units','normalized',...
        'pos',[sx 0.8 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
    annotation('textbox',[sx 0.77 0.2 hei],'String','MM with eff=0.95','EdgeColor','none','FontSize',fontsizesmall,'FontName',fontname);
    %
    uicontrol('Style','Radio','String','Multiple regression (2)','Units','normalized',...
        'pos',[sx 0.7 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
    annotation('textbox',[sx 0.67 0.2 hei],'String','MM with eff=0.90','EdgeColor','none','FontSize',fontsizesmall,'FontName',fontname);
    %
    uicontrol('Style','Radio','String','Hawkins','Units','normalized',...
        'pos',[sx 0.6 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
    annotation('textbox',[sx 0.57 0.2 hei],'String','S with bdp=0.5','EdgeColor','none','FontSize',fontsizesmall,'FontName',fontname);
    %
    uicontrol('Style','Radio','String','Stack loss (y)','Units','normalized',...
        'pos',[sx 0.5 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
    annotation('textbox',[sx 0.47 0.2 hei],'String','LTS','EdgeColor','none','FontSize',fontsizesmall,'FontName',fontname);
    %
    uicontrol('Style','Radio','String','Stack loss (sqrt y)','Units','normalized',...
        'pos',[sx 0.4 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
    annotation('textbox',[sx 0.37 0.2 hei],'String','LTS','EdgeColor','none','FontSize',fontsizesmall,'FontName',fontname);
    %
    str = sprintf(['\n'...
        'To explore how dynamic brushing works, select one of the radio buttons on the left. \n\n' ...
        'The index plot of robust residuals will appear automatically. \n\n'...
        'The units associated with the residuals you select in the index plot of robust residuals will also be automatically highlighted in the yXplot. \n\n'...
        'To produce an index plot of robust residuals with your own data, type:\n\n\n' ...
        '        [out]=LXS(y,X); Least median of squares\n' ...
        '        [out]=LXS(y,X,''lms'',0); Least trimmed squares\n' ...
        '        [out]=Sreg(y,X); S estimator\n' ...
        '        [out]=MMreg(y,X); MM estimator\n\n' ...
        '        resindexplot(out); index plot of robust residuals \n\n\n'...
        'See the help files of the above functions for more information \n\n'...
        'or see file']);
    
    annotation('textbox',[0.3 0 0.7 1],'String',str,'Interpreter','none','FontSize',fontsizebig,'FontName',fontname,'BackgroundColor',backgroundcolor);
    
    import com.mathworks.mlwidgets.html.HTMLRenderer;
    
    % create component
    r = HTMLRenderer;
    % set the text to display
    r.setHtmlText(['<html> <a href="matlab:opentoline(''',examp1,''',5)">example_regression.m</a></html>']);
    % make sure the component is opaque
    r.setOpaque(true);
    % add the component
    javacomponent(r, [290 25 400 35], gcf);
    %[left, bottom, width, height]
    
    % Include at the bottom of the GUI two toggle buttons named close and
    % closeall
    uicontrol('Style','Pushbutton','Units','normalized', ...
        'Position',[0.01 0.07 0.17 0.05], ...
        'Callback','close', 'String','Close current GUI');
    uicontrol('Style','Pushbutton','Units','normalized', ...
        'Position',[0.01 0.01 0.28 0.05], ...
        'Callback','cabc', 'String','Close all plots but current GUI');
    
    % Initialize some button group properties.
    set(h,'SelectionChangeFcn',@brushROBex);
    set(h,'SelectedObject',[]);  % No selection
    set(h,'Visible','on'); % Make the object visible
    
    % Add the logo to the GUI
    imdata = imread('logo.png','BackgroundColor',(240/255)*[1 1 1]);
    %imdata = imread('logo.png');
    %hA=axes('Position',[0.05 0.15 0.2 0.1]);
    
    hA=axes('Position',[0.05 0.15 0.2 0.2],'Layer','top');
    %imdisp(imdata)
    
    image(imdata,'Parent',hA);
    set(hA,'DataAspectRatio',[1 1 1]);
    xlim(hA,[1 64]);
    ylim(hA,[0 65]);
    
    axis('off');
else
    brushROBex([],eventdata)
end

    function brushROBex(~,eventdata)
        %% Examples inside GUI to show how interactive brushing works
        
        if eventdata==1 || (isa(eventdata,'matlab.ui.eventdata.SelectionChangedData') && strcmp(get(eventdata.NewValue,'String'),'Forbes dataset'));
            forbes=load('forbes.txt');
            y=forbes(:,2);
            X=forbes(:,1);
            [out]=LXS(y,X,'nsamp',0,'yxsave',1);
            
            if isa(eventdata,'matlab.ui.eventdata.SelectionChangedData')
                databrush=struct;
                databrush.selectionmode='Rect';
                databrush.persist='';
                databrush.Label='off';
                databrush.RemoveLabels='off';
            else
                databrush='';
            end
            resindexplot(out,'databrush',databrush);
            
            
        elseif eventdata==2 || (isa(eventdata,'matlab.ui.eventdata.SelectionChangedData') && strcmp(get(eventdata.NewValue,'String'),'Multiple regression'))
            multiple_regression=load('multiple_regression.txt');
            y=multiple_regression(:,4);
            X=multiple_regression(:,1:3);
            % LMS using 1000 subsamples
            [out]=MMreg(y,X,'Snsamp',500,'eff',0.95,'yxsave',1);
            
            if isa(eventdata,'matlab.ui.eventdata.SelectionChangedData')
                % plot residuals using brushing
                databrush=struct;
                databrush.bivarfit='';
                databrush.selectionmode='Rect';
                databrush.persist='';
                databrush.Label='off';
                databrush.RemoveLabels='off';
            else
                databrush='';
            end
            
            resindexplot(out,'databrush',databrush);
            
        elseif eventdata==3 || (isa(eventdata,'matlab.ui.eventdata.SelectionChangedData') && strcmp(get(eventdata.NewValue,'String'),'Multiple regression (2)'))
            multiple_regression=load('multiple_regression.txt');
            y=multiple_regression(:,4);
            X=multiple_regression(:,1:3);
            % LMS using 1000 subsamples
            [out]=MMreg(y,X,'Snsamp',500,'eff',0.90,'yxsave',1);
            
            if isa(eventdata,'matlab.ui.eventdata.SelectionChangedData')
                % plot residuals using brushing
                databrush=struct;
                databrush.bivarfit='';
                databrush.selectionmode='Rect';
                databrush.persist='';
                databrush.Label='off';
                databrush.RemoveLabels='off';
            else
                databrush='';
            end
            
            resindexplot(out,'databrush',databrush,'numlab',{6});
            
        elseif eventdata==4 || (isa(eventdata,'matlab.ui.eventdata.SelectionChangedData') && strcmp(get(eventdata.NewValue,'String'),'Hawkins'))
            %% Hawkins data
            hawkins=load('hawkins.txt');
            y=hawkins(:,9);
            X=hawkins(:,1:8);
            [out]=Sreg(y,X,'nsamp',500,'yxsave',1);
            
            if isa(eventdata,'matlab.ui.eventdata.SelectionChangedData')
                databrush=struct;
                databrush.bivarfit='';
                databrush.selectionmode='Brush';
                databrush.persist='';
                databrush.Label='off';
                databrush.RemoveLabels='off';
            else
                databrush='';
            end
            
            resindexplot(out,'databrush',databrush);
            
        elseif eventdata==5 || (isa(eventdata,'matlab.ui.eventdata.SelectionChangedData') && strcmp(get(eventdata.NewValue,'String'),'Stack loss (y)'))
            %% Stack loss data (original scale)
            stack_loss=load('stack_loss.txt');
            y=stack_loss(:,4);
            X=stack_loss(:,1:3);
            [out]=LXS(y,X,'nsamp',0,'lms',0,'yxsave',1);
            
            if isa(eventdata,'matlab.ui.eventdata.SelectionChangedData')
                databrush=struct;
                databrush.bivarfit='2';
                databrush.selectionmode='Rect';
                databrush.persist='';
                databrush.Label='off';
                databrush.RemoveLabels='off';
            else
                databrush='';
            end
            
            resindexplot(out,'databrush',databrush);
            
        elseif eventdata==6 || (isa(eventdata,'matlab.ui.eventdata.SelectionChangedData') && strcmp(get(eventdata.NewValue,'String'),'Stack loss (sqrt y)'))
            %% Stack loss data (sqrt scale)
            stack_loss=load('stack_loss.txt');
            y=sqrt(stack_loss(:,4));
            X=stack_loss(:,1:3);
            [out]=LXS(y,X,'nsamp',0,'lms',0,'yxsave',1);
            
            if isa(eventdata,'matlab.ui.eventdata.SelectionChangedData')
                databrush=struct;
                databrush.bivarfit='2';
                databrush.selectionmode='Rect';
                databrush.persist='';
            else
                databrush='';
            end
            
            resindexplot(out,'databrush',databrush);
            
        else
        end
        
        %%         % Now position the plots in particular areas of the screen
        scrsz = get(0,'ScreenSize');
        % Check which plots are open
        % ScreenSize is a four-element vector: [left, bottom, width, height]:
        width=scrsz(3)/4;
        
        % Check if figure containing mdr is present and get its handle
        plmdr=((findobj('type','figure','Tag','pl_mdr')));
        
        % Check if figure containing residuals is present and get its handle
        plresfwd=((findobj('type','figure','Tag','pl_resfwd')));
        
        % Check if yX plot is present and get its handle
        plyX=((findobj('type','figure','Tag','pl_yX')));
        
        if ~isempty(plmdr) && ~isempty(plresfwd) && ~isempty(plyX)
            % brushing resindexplot generates a yXplot which is
            % automatically positioned. However, if the mdr plot or the
            % forward plot of residuals are presnt, position of the plots
            % is set as follows.
            set(plresfwd,'Position',[10 scrsz(4)/10  width scrsz(4)/3])
            set(plmdr,'Position',[(width+20) scrsz(4)/10  width scrsz(4)/3])
            set(plyX,'Position',[(2*width+40) scrsz(4)/10  width scrsz(4)/3])
        end
        
    end


a=version;
if str2double(a(1))>=8
    stri='Detailed information about the datasets used in this GUI can be found in the <a href="matlab: doc -classic">section USER GUIDE of the FSDA html help system </a>';
else
    stri='Detailed information about the datasets used in this GUI can be found <a href="matlab: docsearchFS(''datasets_reg'')">here</a>';
end
disp(stri)

end