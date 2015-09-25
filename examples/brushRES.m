function brushRES
%brushRES displays a GUI where it is possible to brush steps from the monitoring residuals plot 
%and to see the corresponding units highlighted in other plots
%
%<a href="matlab: docsearchFS('brushRES')">Link to the help page for this function</a>
%
% See also: brushFAN, brushROB
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('brushRES')">Link to the help page for this function</a>
% Last modified 06-Feb-2015

% other open demos create problems. Delete them before starting this new one. 
delete(((findobj('type','figure','Tag','demo'))));

figure; % Create the button group.

h = uibuttongroup('visible','off','Position',[0 0 .3 1]);

% Set name of the GUI
set(gcf,'Name', 'Example of robust dynamic brushing starting from the monitoring residuals plot', 'NumberTitle', 'off','Tag','demo');

wid=0.95;
sx=0.05;
hei=0.04;
fontsizebig = 10;
fontsizesmall = 9;
fontname = 'Lucida';
backgroundcolor = 'w';

% Include on the left the radio buttons
uicontrol('Style','Radio','String','Forbes dataset','Units','normalized',...
    'pos',[sx 0.9 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
uicontrol('Style','Radio','String','Multiple regression','Units','normalized',...
    'pos',[sx 0.8 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
uicontrol('Style','Radio','String','Multiple regression 2','Units','normalized',...
    'pos',[sx 0.7 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
uicontrol('Style','Radio','String','Hawkins','Units','normalized',...
    'pos',[sx 0.6 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
uicontrol('Style','Radio','String','Stack loss (y)','Units','normalized',...
    'pos',[sx 0.5 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);
uicontrol('Style','Radio','String','Stack loss (sqrt y)','Units','normalized',...
    'pos',[sx 0.4 wid hei],'parent',h,'HandleVisibility','off','FontSize',fontsizesmall,'FontName',fontname);

examp=which('examples_regression.m');
examp1=strrep(examp,'\','\\');

str = sprintf(['\n'...
    'To explore how dynamic brushing works, select one of the radio buttons on the left. \n\n' ...
    'The monitoring of the residuals plot together with the minimum deletion residual plot will appear automatically. \n\n'...
    'The points associated with the trajectories you select in the plot of scaled residuals will also be automatically highlighted in the plots of minimum deletion residual and in the scatter plot matrix. \n\n'...
    'To monitor the residuals and the minimum deletion residual with your own data, type:\n\n\n' ...
    '        [out]=LXS(y,X); to find a robust initial subset\n' ...
    '        [out]=FSReda(y,X,out.bs);  To run the forward search \n'...
    '        mdrplot(out); To plot the minimum deletion residual  \n'...
    '        resfwdplot(out); To plot the scaled squared residuals \n\n\n'...
    'See the help files of the above functions for more information \n\n'...
    'or see file']);


annotation('textbox',[0.3 0 0.7 1],'String',str,'FontSize',fontsizebig,'FontName',fontname,'BackgroundColor',backgroundcolor,'Interpreter','none');

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
set(h,'SelectionChangeFcn',@brushRESex);
set(h,'SelectedObject',[]);  % No selection
set(h,'Visible','on'); % Make the object visible

% Add the logo to the GUI
imdata = imread('logo.png','BackgroundColor',(240/255)*[1 1 1]);

hA=axes('Position',[0.05 0.15 0.2 0.2],'Layer','top');

image(imdata,'Parent',hA);
set(hA,'DataAspectRatio',[1 1 1]);
xlim(hA,[1 64]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(hA,[0 65]);

axis('off')

    function brushRESex(~,eventdata)
        %% Examples inside GUI to show how interactive brushing works
        
        % Set graphical options
        fsiztitl=12; % Font size of title and of the label of the axes
        SizeAxesNum=12; % Font size of the numbers on the axes
        
        if strcmp(get(eventdata.NewValue,'String'),'Forbes dataset');
            forbes=load('forbes.txt');
            y=forbes(:,2);
            X=forbes(:,1);
            [out]=LXS(y,X,'nsamp',0);
            [out]=FSReda(y,X,out.bs);
            
            % Plot minimum deletion residual
            mdrplot(out,'xlimx',[6 17],'ylimy',[0 13],'FontSize',fsiztitl,'SizeAxesNum',SizeAxesNum);
            
            % Plot monitoring of scaled residuals
            databrush=struct;
            databrush.bivarfit='i1';
            databrush.selectionmode='Rect';
            databrush.persist='';
            databrush.Label='on';
            databrush.RemoveLabels='off';
            standard=struct;
            standard.SizeAxesNum=SizeAxesNum;
            cascade;
            resfwdplot(out,'databrush',databrush,'standard',standard);
            
        elseif strcmp(get(eventdata.NewValue,'String'),'Multiple regression')
            multiple_regression=load('multiple_regression.txt');
            y=multiple_regression(:,4);
            X=multiple_regression(:,1:3);
            % LMS using 1000 subsamples
            [out]=LXS(y,X,'nsamp',1000);
            % Forward Search
            [out]=FSReda(y,X,out.bs);
            out1=out;
            % Create scaled squared residuals
            out1.RES=out.RES.^2;
            % plot minimum deletion residual with personalized options
            mdrplot(out,'ylimy',[1 4.2],'xlimx',[10 60],'FontSize',fsiztitl,'SizeAxesNum',SizeAxesNum,'lwdenv',2);
            
            % plot the scaled residuals using brushing
            databrush=struct;
            databrush.bivarfit='';
            databrush.selectionmode='Rect';
            databrush.persist='';
            databrush.Label='on';
            databrush.RemoveLabels='off';
            standard=struct;
            standard.SizeAxesNum=SizeAxesNum;
            cascade;
            resfwdplot(out,'databrush',databrush,'standard',standard);
            
        elseif strcmp(get(eventdata.NewValue,'String'),'Multiple regression 2')
            multiple_regression=load('multiple_regression.txt');
            y=multiple_regression(:,4);
            X=multiple_regression(:,1:3);
            % LMS using 1000 subsamples
            [out]=LXS(y,X,'nsamp',10000);
            % Forward Search
            [out]=FSReda(y,X,out.bs);
            out1=out;
            % Create scaled squared residuals
            out1.RES=out.RES.^2;
            
            seltyp={'--' '-' '-.' ':'};
            selcolor={'b';'g';'c';'m';'y';'k'}; % Specify the colors for the lines in
            standard=struct;
            standard.SizeAxesNum=SizeAxesNum;
            standard.Color=selcolor;
            standard.xlim=[12 60];
            standard.ylim=[0 30];
            standard.LineStyle=seltyp;
            
            databrush=struct;
            databrush.bivarfit='i1';
            databrush.selectionmode='Rect';
            databrush.persist='';
            databrush.Label='on';
            databrush.RemoveLabels='on';
            
            fground=struct;
            fground.fthresh=15;
            fground.LineWidth=3;
            cascade;
            resfwdplot(out1,'databrush',databrush,'standard',standard,'fground',fground);
            
        elseif strcmp(get(eventdata.NewValue,'String'),'Hawkins')
            %% Hawkins data
            hawkins=load('hawkins.txt');
            y=hawkins(:,9);
            X=hawkins(:,1:8);
            [out]=LXS(y,X,'lms',0,'nsamp',2000);
            [out]=FSReda(y,X,out.bs);
            % Plot minimum deletion residual using personalized graphical options
            mdrplot(out,'ylimy',[1 8],'xlimx',[25 128],'FontSize',fsiztitl,'SizeAxesNum',SizeAxesNum,'lwdenv',2);
            out1=out;
            out1.RES=out.RES;
            selcolor={'b'};
            seltyp={'--' '-.' ':'};
            
            standard=struct;
            standard.SizeAxesNum=SizeAxesNum;
            standard.Color=selcolor;
            standard.xlim=[21 128];
            standard.ylim=[-4.5 4.5];
            standard.LineStyle=seltyp;
            
            databrush=struct;
            databrush.bivarfit='';
            databrush.selectionmode='Brush';
            databrush.persist='';
            databrush.Label='off';
            databrush.RemoveLabels='off';
            
            fground=struct;
            fground.fthresh=1.8;
            fground.LineWidth=3;
            fground.LineStyle={'-'};
            
            cascade;
            resfwdplot(out1,'databrush',databrush,'standard',standard,'fground',fground);
            
        elseif strcmp(get(eventdata.NewValue,'String'),'Stack loss (y)')
            %% Stack loss data (original scale)
            stack_loss=load('stack_loss.txt');
            y=stack_loss(:,4);
            X=stack_loss(:,1:3);
            [out]=LXS(y,X,'nsamp',0);
            [out]=FSReda(y,X,out.bs,'init',5);
            mdrplot(out,'ylimy',[0.5 5],'xlimx',[5 21]);
            
            databrush=struct;
            databrush.bivarfit='2';
            databrush.selectionmode='Rect';
            databrush.persist='';
            databrush.Label='on';
            databrush.RemoveLabels='off';
            cascade;
            resfwdplot(out,'databrush',databrush);
            
        elseif strcmp(get(eventdata.NewValue,'String'),'Stack loss (sqrt y)')
            %% Stack loss data (sqrt scale)
            stack_loss=load('stack_loss.txt');
            y=sqrt(stack_loss(:,4));
            X=stack_loss(:,1:3);
            [out]=LXS(y,X,'nsamp',0);
            [out]=FSReda(y,X,out.bs,'init',5);
            mdrplot(out,'ylimy',[0.5 5],'xlimx',[5 21]);
            databrush=struct;
            databrush.bivarfit='2';
            databrush.selectionmode='Rect';
            databrush.persist='on';
            databrush.Label='on';
            databrush.RemoveLabels='off';
            cascade;
            resfwdplot(out,'databrush',databrush);
            
        else
        end
        
        %% Now position the plots in particular areas of the screen
        
        % Ensure root units are pixels and get the size of the screen

        set(0,'Units','pixels') 
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
        
        % brushing forward plot of residuals generates a yXplot which is
        % automatically positioned. However, if the mdr plot is also are
        % present, position of the plots is set as follows.
        if ~isempty(plmdr) && ~isempty(plresfwd) && ~isempty(plyX)
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

% OLD to delete
%    'See file <a href="matlab: opentoline(' examp1 ',27)">example_regression.m</a>  or the help files of the above functions for more information']);
%    'See file "examples_regression.m" or the help files of the above functions for more information']);
%
%    stri=['See file <a href="matlab: opentoline(' examp1 ',27)">example_regression.m</a>  or the help files of the above functions for more information'];
%    sprintf(stri)
% annotation('textbox',[0.3 0 0.7 1],'String',str,'Interpreter','none','FontSize',fontsizebig,'FontName',fontname,'BackgroundColor',backgroundcolor);
%stri=['See file <a href="matlab: opentoline(' examp1 ',27)">example_regression.m</a>  or the help files of the above functions for more information'];
%disp(stri)
%
%
%<a href="matlab: docsearchFS('brushRES')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    % Interactive_example
    % Display a GUI where it is possible to brush steps from the monitoring
    % residuals plot and to see the corresponding units highlighted in
    % other plots
    brushRES
%}